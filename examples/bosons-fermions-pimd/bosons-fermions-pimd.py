"""
Path integral molecular dynamics of bosons and fermions
========================================================

:Authors: Barak Hirshberg `@bhirshberg <https://github.com/bhirshberg>`_

This example shows how to run path integral molecular dynamics (PIMD)
simulations of **indistinguishable** particles -- bosons and fermions -- with
``i-PI``, and how to analyze the output. It is a modernised version of the
"Bosonic and Fermionic PIMD" module of the `PIQM 2023 tutorial
<https://github.com/i-pi/piqm2023-tutorial>`_, updated to i-PI 3.x.

We simulate three non-interacting particles in a three-dimensional harmonic
trap (:math:`\\hbar\\omega_0 = 3\\,\\mathrm{meV}`) and compute the average energy
for four kinds of statistics:

* three distinguishable particles (the baseline),
* three bosons,
* a mixture of two bosons and one distinguishable particle,
* three fermions.

The forces come from the built-in ``harmonic`` potential of i-PI's Python
driver, so the whole recipe installs from PyPI -- no compiled driver required.
"""

import concurrent.futures
import os
import re
import shutil
import subprocess
import tempfile
import time

import matplotlib.pyplot as plt
import numpy as np

import analysis


# %%
# Indistinguishable particles and ring-polymer exchange
# -----------------------------------------------------
#
# In path integral molecular dynamics each quantum particle is represented by a
# ring polymer of :math:`P` beads. For **distinguishable** particles every ring
# polymer is closed on itself. Quantum **exchange** between identical particles
# is captured by also allowing ring polymers of different particles to be
# connected into longer rings. Averaging over all such connectivities with the
# appropriate weights yields bosonic statistics; fermions additionally require a
# sign :math:`(-1)^{\\ldots}` for every pair exchange.
#
# In i-PI this is switched on with a single tag inside ``<normal_modes>``:
#
# .. code-block:: xml
#
#     <normal_modes propagator='bab'>
#         <nmts> 10 </nmts>
#         <bosons> [0, 1, 2] </bosons>
#     </normal_modes>
#
# ``<bosons>`` holds the (zero-based) indices of the atoms that participate in
# bosonic exchange. An empty list recovers distinguishable particles; listing a
# subset gives a mixture. The exchange spring potential is evaluated with the
# quadratic-scaling algorithm of Hirshberg *et al.*
# (`PNAS 2019 <https://doi.org/10.1073/pnas.1913365116>`_) and Feldman &
# Hirshberg (`JCP 2023 <https://doi.org/10.1063/5.0173749>`_).


# %%
# Running an i-PI simulation from Python
# --------------------------------------
#
# i-PI uses a client-server model: the ``i-pi`` server evolves the ring polymers
# and a *driver* computes the forces, talking to each other over a socket.
# ``-o`` passes the spring constant
# :math:`k = 1.216\\times10^{-8}\\,\\mathrm{Ha/Bohr^2}` matching
# :math:`\\hbar\\omega_0 = 3\\,\\mathrm{meV}`.
#
# For the forces we use the harmonic driver that ships with i-PI. The
# pure-Python ``i-pi-py_driver`` (mode ``harmonic``) always works and needs no
# compilation; if the faster compiled Fortran driver ``i-pi-driver`` (mode
# ``harm3d``) is on the ``PATH`` we use it instead -- it is a few times faster
# and gives identical results. The helper below picks whichever is available,
# launches both processes, waits, and stores the output of each run under a
# case-specific name (``data_<tag>.out``).

SPRING_CONSTANT = "1.21647924e-8"  # Ha/Bohr^2


def driver_command(address):
    """Return the force-driver command, preferring the compiled f90 driver."""
    if shutil.which("i-pi-driver") is not None:  # compiled Fortran driver
        return ["i-pi-driver", "-m", "harm3d", "-o", SPRING_CONSTANT, "-u", "-a", address]
    return ["i-pi-py_driver", "-m", "harmonic", "-o", SPRING_CONSTANT, "-u", "-a", address]


def run_case(tag, xml, address):
    """Run one i-PI + driver simulation and return the output-file path."""
    # start from a clean slate so i-PI does not refuse to overwrite outputs
    for stale in ["data.out"] + [f for f in os.listdir(".") if f.startswith("data.pos")]:
        if os.path.exists(stale):
            os.remove(stale)
    sock = f"/tmp/ipi_{address}"
    if os.path.exists(sock):
        os.remove(sock)

    ipi = subprocess.Popen(["i-pi", f"data/{xml}"], stdout=subprocess.DEVNULL)
    time.sleep(2)  # give i-PI time to open the socket
    driver = subprocess.Popen(driver_command(address), stdout=subprocess.DEVNULL)
    ipi.wait()
    driver.wait()

    out = f"data_{tag}.out"
    os.replace("data.out", out)
    return out


def _run_seeded(xml, seed, address):
    """Run one trajectory with a given random seed in an isolated temp dir and
    return the path to its ``data.out`` (used for the multi-trajectory demo)."""
    tmp = tempfile.mkdtemp(prefix="ipi_traj_")
    text = open(f"data/{xml}").read()
    text = re.sub(r"<seed>.*?</seed>", f"<seed> {seed} </seed>", text, flags=re.S)
    text = re.sub(r"<address>.*?</address>", f"<address>{address}</address>", text, flags=re.S)
    open(os.path.join(tmp, "input.xml"), "w").write(text)

    sock = f"/tmp/ipi_{address}"
    if os.path.exists(sock):
        os.remove(sock)
    ipi = subprocess.Popen(["i-pi", "input.xml"], cwd=tmp,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(120):  # wait for the socket before starting the driver
        if os.path.exists(sock) or ipi.poll() is not None:
            break
        time.sleep(0.5)
    drv = subprocess.Popen(driver_command(address), cwd=tmp,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    ipi.wait()
    drv.wait()
    return os.path.join(tmp, "data.out")


def run_fermion_ensemble(n_traj, skip_rows):
    """Run ``n_traj`` short fermionic trajectories (different seeds) in parallel
    and combine them with the sign-weighted estimator. Returns
    ``(mean, error, n_eff, mean_sign)``."""
    seeds = [4001 + 137 * i for i in range(n_traj)]
    workers = min(n_traj, max(1, (os.cpu_count() or 2) // 2))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        outs = list(pool.map(
            lambda i: _run_seeded("input_3fermions.xml", seeds[i], f"bf-fmulti-{i}"),
            range(n_traj),
        ))
    E, W, signs = [], [], []
    for out in outs:
        e_j, w_j, s_j = analysis.fermionic_trajectory_estimate(out, skip_rows)
        E.append(e_j)
        W.append(w_j)
        signs.append(s_j)
        shutil.rmtree(os.path.dirname(out), ignore_errors=True)
    mean, error, n_eff = analysis.weighted_average(E, W)
    return mean, error, n_eff, float(np.mean(signs))


# %%
# Distinguishable particles (the baseline)
# ----------------------------------------
#
# We first run three *distinguishable* particles (``<bosons> []``). For a bound
# system we can use two kinetic-energy estimators: the primitive/quantum virial
# (``virial_fq``) and the centroid virial (``kinetic_cv``). Both should agree
# with the analytical result. We discard the first part of the trajectory as
# thermalisation.
#
# .. note::
#    ``total_steps`` is deliberately small (1000--3000 steps, with a large
#    ``timestep`` of 10 fs allowed by the very soft trap) so the whole recipe
#    runs in a few minutes -- the results below are therefore only qualitative.
#    For converged averages you should use at least :math:`10^6` steps; the
#    exact analytical values are printed alongside for reference.

skip = 200  # rows to discard as thermalisation (properties are written every step)

dist_out = run_case("3dist", "input_3dist.xml", "bf-3dist")
dist = analysis.mean_energies(dist_out, skip)
dist_ref = analysis.analytical_energy(17.4, "dist")

print("Distinguishable particles (T = 17.4 K)")
print(f"  total energy, virial estimator   : {dist['total']:.7f} Ha")
print(f"  kinetic, centroid-virial          : {dist['kinetic_cv']:.7f} Ha")
print(f"  analytical total energy           : {dist_ref:.7f} Ha")


# %%
# Bosons
# ------
#
# Turning the three particles into bosons only requires ``<bosons> [0, 1, 2]``.
# Exchange lowers the energy relative to distinguishable particles.
#
# .. warning::
#    The **centroid-virial** kinetic estimator is *not* valid under bosonic
#    exchange, so we do not request it here. The primitive virial (valid for
#    bound systems) and the potential estimator remain correct.

bos_out = run_case("3bosons", "input_3bosons.xml", "bf-3bosons")
bos = analysis.mean_energies(bos_out, skip)
bos_ref = analysis.analytical_energy(17.4, "bosonic")

print("Three bosons (T = 17.4 K)")
print(f"  total energy, virial estimator   : {bos['total']:.7f} Ha")
print(f"  analytical total energy           : {bos_ref:.7f} Ha")


# %%
# A mixture: two bosons and one distinguishable particle
# ------------------------------------------------------
#
# Listing only atoms 0 and 1 in ``<bosons>`` makes them exchange with each
# other while atom 2 stays distinguishable.

mix_out = run_case("2bosons1dist", "input_2bosons1dist.xml", "bf-2bosons1dist")
mix = analysis.mean_energies(mix_out, skip)
mix_ref = analysis.analytical_energy(17.4, "mixed")

print("Two bosons + one distinguishable particle (T = 17.4 K)")
print(f"  total energy, virial estimator   : {mix['total']:.7f} Ha")
print(f"  analytical total energy           : {mix_ref:.7f} Ha")


# %%
# Fermions: reweighting and the sign problem
# ------------------------------------------
#
# Fermionic averages are obtained by simulating **bosons** and reweighting each
# sample by the fermionic sign :math:`s`,
#
# .. math::
#
#    \\langle A \\rangle_F = \\frac{\\langle A\\, s \\rangle}{\\langle s \\rangle},
#
# where the sign follows from the recursive configuration weight
#
# .. math::
#
#    W^{(N)} = \\frac{1}{N}\\sum_{k=1}^{N} \\xi^{k-1}
#              e^{-\\beta E_N^{(k)}} \\, W^{(N-k)}, \\qquad W^{(0)} = 1,
#
# with :math:`\\xi = -1` for fermions. In the 2023 tutorial this weight was
# evaluated by a hand-written module; **i-PI 3.x computes it for us** and writes
# the sign to the ``fermionic_sign`` column, so we just request that property
# and reweight in two lines (see ``analysis.reweighted_fermionic_energy``).
#
# When the average sign approaches zero the estimator becomes hard to converge
# -- the fermionic **sign problem**. This is why we use a higher temperature
# (30 K) and fewer beads for the fermionic run.
#
# .. admonition:: The exact reference value
#
#    The exact three-fermion energy at 30 K is **1.053 mHa**, computed by
#    ``analysis.analytical_energy`` from the canonical partition-function
#    recursion (elementary symmetric polynomial, :math:`\\xi=-1`). An earlier
#    version of this tutorial used an incorrect hard-coded closed form that gave
#    0.912 mHa; the value here is confirmed independently by brute-force
#    enumeration of the three-fermion states. A single short run below is only
#    illustrative -- see the two sections that follow for a properly sampled
#    comparison.

fer_out = run_case("3fermions", "input_3fermions.xml", "bf-3fermions")
mean_sign, fer_energy = analysis.reweighted_fermionic_energy(fer_out, skip)
fer_ref = analysis.analytical_energy(30.0, "fermionic")

print("Three fermions (T = 30 K, single short run)")
print(f"  average sign <s>                  : {mean_sign:.4f}")
print(f"  fermionic total energy            : {fer_energy:.7f} Ha")
print(f"  analytical total energy           : {fer_ref:.7f} Ha")


# %%
# Don't trust one fermionic run -- average several
# ------------------------------------------------
#
# That single fermionic number is almost meaningless: the average sign is small,
# so the reweighted energy has a huge variance and one short run can land far
# from the truth. The honest thing is to run **several independent
# trajectories** (different random seeds) and combine them. Because the sign
# varies between trajectories, we use the sign-*weighted* estimator explained in
# the "Error estimation" section below (not a plain average). Here we run a
# handful of short trajectories in parallel -- still only a couple of minutes on
# a laptop:

fer_mean, fer_err, n_eff, ens_sign = run_fermion_ensemble(n_traj=8, skip_rows=skip)

print("Three fermions (T = 30 K, 8 short trajectories, sign-weighted)")
print(f"  average sign <s>                  : {ens_sign:.3f}")
print(f"  weighted fermionic energy         : {fer_mean * 1e3:.4f} +/- {fer_err * 1e3:.4f} mHa")
print(f"  analytical total energy           : {fer_ref * 1e3:.4f} mHa")
print(f"  (effective sample size n_eff = {n_eff:.1f} of 8)")

# %%
# The error bar is large -- often 5-10% of the value -- but it now *brackets* the
# analytical result. That large-but-honest error bar is the fingerprint of the
# fermionic **sign problem**: the numbers make sense, they are just noisy, and
# tightening them needs far more sampling (see the reference study). Contrast
# this with the distinguishable/bosonic runs, whose single short runs are already
# close because there is no sign to fight.


# %%
# Summary: exchange orders the energies
# -------------------------------------
#
# Bosonic exchange *lowers* the energy while fermionic antisymmetry *raises* it,
# so at fixed temperature :math:`E_\\mathrm{bosons} < E_\\mathrm{dist} <
# E_\\mathrm{fermions}`. The distinguishable/bosonic bars are single short runs
# (already close to the exact value); the fermionic bar is the 8-trajectory
# sign-weighted estimate with its error bar.

labels = ["bosons", "dist", "2B+1D", "fermions*"]
sim = [bos["total"], dist["total"], mix["total"], fer_mean]
sim_err = [0.0, 0.0, 0.0, fer_err]
ref = [bos_ref, dist_ref, mix_ref, fer_ref]

x = np.arange(len(labels))
fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
ax.bar(x - 0.2, np.array(sim) * 1e3, 0.4, yerr=np.array(sim_err) * 1e3,
       capsize=3, label="PIMD")
ax.bar(x + 0.2, np.array(ref) * 1e3, 0.4, label="exact")
ax.set_xticks(x, labels)
ax.set_ylabel("total energy / mHa")
ax.set_title("Three particles in a harmonic trap")
ax.legend()
ax.text(0.02, 0.02, "*fermions at 30 K (8 traj), others at 17.4 K (1 run)",
        transform=ax.transAxes, fontsize=8, color="gray")


# %%
# Error estimation for fermions: weight trajectories by the sign
# --------------------------------------------------------------
#
# For fermions a *single* run gives one number with no error bar, and naively
# averaging the per-trajectory reweighted energies
# :math:`E_j = \\langle \\varepsilon\\, s\\rangle_j / \\langle s\\rangle_j` with a
# plain :math:`\\mathrm{std}/\\sqrt{M}` is **biased**: a trajectory that happened
# to sample mostly low-weight (small-sign) configurations yields a wild ratio,
# yet counts as much as a well-converged one.
#
# The fix (SI of Hirshberg, Invernizzi & Parrinello, *J. Chem. Phys.* **152**,
# 171102 (2020)) is to weight each trajectory by its total weight in the
# reweighted ensemble, :math:`W_j = \\sum_\\mathrm{steps} s` (the sum of the
# instantaneous signs), and to use an *effective* sample size:
#
# .. math::
#
#    \\bar E_F = \\frac{\\sum_j W_j E_j}{\\sum_j W_j}, \\qquad
#    n_\\mathrm{eff} = \\frac{\\left(\\sum_j W_j\\right)^2}{\\sum_j W_j^2}, \\qquad
#    \\sigma_E^2 = \\frac{n_\\mathrm{eff}}{n_\\mathrm{eff}-1}
#                 \\frac{\\sum_j W_j (E_j - \\bar E_F)^2}{\\sum_j W_j},
#
# with statistical error :math:`\\sigma_E/\\sqrt{n_\\mathrm{eff}}`. When all
# weights are equal this reduces to the ordinary mean and
# :math:`\\mathrm{std}/\\sqrt{M}`, which is why the bosonic rows above need no
# special treatment. Both estimators are implemented in ``analysis.py``
# (:func:`weighted_average`, :func:`fermionic_trajectory_estimate`) -- the same
# functions the 8-trajectory demo above used.
#
# Running many more trajectories (here 20 of 5000 steps) makes the difference
# between the two estimators obvious (exact value 1.053 mHa):
#
# ==========================  =====================  =============================
# estimator                   mean (mHa)             error (mHa)
# ==========================  =====================  =============================
# naive mean-of-ratios        1.45  (biased high)    :math:`\\pm 0.32`
# SI weighted, assuming n=M    1.11                   :math:`\\pm 0.078`
# SI weighted, using n_eff     1.11                   :math:`\\pm 0.087`  (n_eff=16 of 20)
# ==========================  =====================  =============================
#
# Two lessons: the weighting **removes the bias** in the mean (the naive
# mean-of-ratios sits at 1.45, the weighted estimate at 1.11), and using the
# effective sample size gives an **honest, ~12% larger** error bar
# (:math:`\\sqrt{M/n_\\mathrm{eff}}`) than pretending all :math:`M` trajectories
# are equally informative. The weighted estimate 1.11 :math:`\\pm` 0.09 agrees
# with the exact 1.053 mHa within its error bar -- there is *no* mysterious
# fermionic discrepancy once (i) the correct benchmark is used and (ii) the
# statistical error is estimated properly.
#
# The runs in this tutorial are deliberately short and use only 12 beads for the
# fermions -- enough to "make sense", not to be tightly converged. For
# production you would simply run more (and longer) trajectories and, if needed,
# increase the number of beads; the estimators above are exactly what you would
# use.


# %%
# Assignment: energy vs temperature
# ---------------------------------
#
# A classic exercise is to plot the mean total energy against
# :math:`\\beta\\hbar\\omega_0` for distinguishable particles, bosons and
# fermions. The analytical curves (cheap, no simulation needed) look like this;
# reproducing the points with converged PIMD runs is left as an exercise.

omega = analysis.omega0()
bhw = np.linspace(0.5, 6.0, 60)
temps = omega / (bhw * analysis.KELVIN_TO_HARTREE)  # T for each beta*hbar*omega

fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
for sys_type, style in [("dist", "-"), ("bosonic", "--"), ("fermionic", ":")]:
    e = [analysis.analytical_energy(T, sys_type) / omega for T in temps]
    ax.plot(bhw, e, style, label=sys_type)
ax.set_xlabel(r"$\beta \hbar \omega_0$")
ax.set_ylabel(r"$\langle E \rangle \,/\, \hbar\omega_0$")
ax.set_title("Analytical energies, three particles")
ax.legend()

# %%
#
# **References**
#
# * B. Hirshberg, V. Rizzi, M. Parrinello,
#   `PNAS 116, 21445 (2019) <https://doi.org/10.1073/pnas.1913365116>`_
# * Y. M. Y. Feldman, B. Hirshberg,
#   `J. Chem. Phys. 159, 154107 (2023) <https://doi.org/10.1063/5.0173749>`_
# * Original PIQM 2023 tutorial:
#   https://github.com/i-pi/piqm2023-tutorial
