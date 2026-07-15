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

import os
import subprocess
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
# and a *driver* computes the forces, talking to each other over a socket. Here
# the driver is ``i-pi-py_driver`` with the ``harmonic`` PES; ``-o`` passes the
# spring constant :math:`k = 1.216\\times10^{-8}\\,\\mathrm{Ha/Bohr^2}` matching
# :math:`\\hbar\\omega_0 = 3\\,\\mathrm{meV}`.
#
# The helper below launches both processes, waits for them to finish, and stores
# the output of each run under a case-specific name (``data_<tag>.out``).

SPRING_CONSTANT = "1.21647924e-8"  # Ha/Bohr^2


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
    driver = subprocess.Popen(
        ["i-pi-py_driver", "-m", "harmonic", "-o", SPRING_CONSTANT,
         "-u", "-a", address],
        stdout=subprocess.DEVNULL,
    )
    ipi.wait()
    driver.wait()

    out = f"data_{tag}.out"
    os.replace("data.out", out)
    return out


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

fer_out = run_case("3fermions", "input_3fermions.xml", "bf-3fermions")
mean_sign, fer_energy = analysis.reweighted_fermionic_energy(fer_out, skip)
fer_ref = analysis.analytical_energy(30.0, "fermionic")

print("Three fermions (T = 30 K)")
print(f"  average sign <s>                  : {mean_sign:.4f}")
print(f"  fermionic total energy            : {fer_energy:.7f} Ha")
print(f"  analytical total energy           : {fer_ref:.7f} Ha")


# %%
# Summary: exchange orders the energies
# -------------------------------------
#
# Bosonic exchange *lowers* the energy while fermionic antisymmetry *raises* it,
# so at fixed temperature :math:`E_\\mathrm{bosons} < E_\\mathrm{dist} <
# E_\\mathrm{fermions}`. (The short runs above are noisy; the ordering and
# agreement become clean once the simulations are converged.)

labels = ["bosons", "dist", "2B+1D", "fermions*"]
sim = [bos["total"], dist["total"], mix["total"], fer_energy]
ref = [bos_ref, dist_ref, mix_ref, fer_ref]

x = np.arange(len(labels))
fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
ax.bar(x - 0.2, np.array(sim) * 1e3, 0.4, label="PIMD (short run)")
ax.bar(x + 0.2, np.array(ref) * 1e3, 0.4, label="analytical")
ax.set_xticks(x, labels)
ax.set_ylabel("total energy / mHa")
ax.set_title("Three particles in a harmonic trap")
ax.legend()
ax.text(0.02, 0.02, "*fermions at 30 K, others at 17.4 K",
        transform=ax.transAxes, fontsize=8, color="gray")


# %%
# Convergence: averaging over independent trajectories
# ----------------------------------------------------
#
# A single short run is noisy, but the dynamics are unbiased (NVT with a proper
# thermostat), so the average over several **independent** trajectories -- each
# started from a different random seed -- converges towards the exact value, and
# the spread across trajectories gives a proper standard error of the mean.
#
# The script ``reference/run_convergence.py`` runs this study (10 trajectories
# per case, parallelised across CPU cores, using the *same* short settings as
# above). Running it reproduces the following table (energies in mHa):
#
# ===================  ==========================  =================
# case                 trajectory average           analytical
# ===================  ==========================  =================
# 3 distinguishable    0.643 :math:`\\pm` 0.017      0.651
# 3 bosons             0.569 :math:`\\pm` 0.017      0.580
# 2 bosons + 1 dist    0.606 :math:`\\pm` 0.018      0.624
# 3 fermions           see next section             0.912
# ===================  ==========================  =================
#
# (10 trajectories per case; the fermionic case needs the weighted estimator of
# the next section.)
#
# The three distinguishable/bosonic cases now agree with the analytical result
# within about one standard error -- for equal-weight data (bosons, and the
# average sign) averaging over trajectories simply removes the noise of a single
# short run. The **fermionic** row is different and needs more care, as we now
# explain.


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
# (:func:`weighted_average`, :func:`fermionic_trajectory_estimate`) and used by
# ``reference/run_convergence.py``.
#
# For 20 trajectories of 5000 steps the two estimators give (analytical
# 0.912 mHa):
#
# ==========================  =====================  =============================
# estimator                   mean (mHa)             error (mHa)
# ==========================  =====================  =============================
# naive mean-of-ratios        1.45  (biased high)    :math:`\\pm 0.32`
# SI weighted, assuming n=M    1.11                   :math:`\\pm 0.078`
# SI weighted, using n_eff     1.11                   :math:`\\pm 0.087`  (n_eff=16 of 20)
# ==========================  =====================  =============================
#
# Two lessons: the weighting **removes the bias** in the mean (1.45 -> 1.11),
# and using the effective sample size gives an **honest, ~12% larger** error bar
# (:math:`\\sqrt{M/n_\\mathrm{eff}}`) than pretending all :math:`M` trajectories
# are equally informative. The residual gap to the exact 0.912 mHa is physical
# (finite beads and the slow convergence of the sign) -- closing it needs far
# more sampling than a couple-minute tutorial allows, which is the whole point
# of the sign problem.


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
