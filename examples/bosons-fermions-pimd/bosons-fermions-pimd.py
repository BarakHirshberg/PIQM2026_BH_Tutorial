"""
Path integral molecular dynamics of bosons and fermions
========================================================

:Authors: Barak Hirshberg `@bhirshberg <https://github.com/bhirshberg>`_

This example shows how to run path integral molecular dynamics (PIMD)
simulations of **indistinguishable** particles -- bosons and fermions -- with
``i-PI``, and how to analyze the output. It is a modernised version of the
"Bosonic and Fermionic PIMD" module of the `PIQM 2023 tutorial
<https://github.com/i-pi/piqm2023-tutorial>`_, updated to i-PI 3.x.

We simulate three non-interacting particles (mass 1) in a three-dimensional
isotropic harmonic trap with :math:`\\hbar\\omega_0 = 3\\,\\mathrm{meV}` (a very
soft trap) and study how quantum statistics changes the average energy. We work
in the natural dimensionless inverse temperature
:math:`\\beta\\hbar\\omega_0 = \\hbar\\omega_0 / k_\\mathrm{B}T`, so the whole
problem is controlled by a single number.

The tutorial is built in three steps:

#. **Bosons.** Switch on bosonic exchange with a one-line tag and trace the
   energy as a function of temperature, comparing against the exact result.
#. **Switching statistics.** See how trivially the *same* input flips between
   distinguishable particles, bosons and mixtures, and compare their energies at
   one temperature.
#. **Fermions.** Obtain fermionic averages by reweighting, meet the **sign
   problem**, and learn how to get an honest error bar by averaging several
   trajectories.

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
# sign :math:`(-1)^{\ldots}` for every pair exchange.
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
# and a *driver* computes the forces, talking to each other over a socket. The
# driver is passed the spring constant
# :math:`k = 1.216\times10^{-8}\,\mathrm{Ha/Bohr^2}`, which corresponds to
# :math:`\hbar\omega_0 = 3\,\mathrm{meV}`.
#
# For the forces we use the harmonic driver that ships with i-PI. The
# pure-Python ``i-pi-py_driver`` (mode ``harmonic``) always works and needs no
# compilation; if the faster compiled Fortran driver ``i-pi-driver`` (mode
# ``harm3d``) is on the ``PATH`` we use it instead -- it is a few times faster
# and gives identical results.
#
# The helper below takes one of the four template inputs in ``data/`` and runs
# it in an isolated temporary directory, optionally overriding the temperature,
# the number of beads, the random seed or the number of steps. This lets us
# drive *every* simulation in the tutorial from a single function.

SPRING_CONSTANT = "1.21647924e-8"  # Ha/Bohr^2, matches hbar*omega0 = 3 meV


def driver_command(address):
    """Return the force-driver command, preferring the compiled f90 driver."""
    if shutil.which("i-pi-driver") is not None:  # compiled Fortran driver
        return [
            "i-pi-driver",
            "-m",
            "harm3d",
            "-o",
            SPRING_CONSTANT,
            "-u",
            "-a",
            address,
        ]
    return [
        "i-pi-py_driver",
        "-m",
        "harmonic",
        "-o",
        SPRING_CONSTANT,
        "-u",
        "-a",
        address,
    ]


def _prepare_input(
    src_xml, address, temp=None, nbeads=None, seed=None, total_steps=None
):
    """Read a template input from ``data/`` and apply the requested overrides."""
    text = open(f"data/{src_xml}").read()
    text = re.sub(
        r"<address>.*?</address>", f"<address>{address}</address>", text, flags=re.S
    )
    if temp is not None:
        # update both the ensemble temperature and the thermal velocity init
        text = re.sub(
            r"(<temperature[^>]*>)\s*[\d.eE+-]+\s*(</temperature>)",
            rf"\g<1> {temp} \g<2>",
            text,
        )
        text = re.sub(
            r"(<velocities[^>]*>)\s*[\d.eE+-]+\s*(</velocities>)",
            rf"\g<1> {temp} \g<2>",
            text,
        )
    if nbeads is not None:
        text = re.sub(r'nbeads="\d+"', f'nbeads="{nbeads}"', text)
    if seed is not None:
        text = re.sub(r"<seed>.*?</seed>", f"<seed> {seed} </seed>", text, flags=re.S)
    if total_steps is not None:
        text = re.sub(
            r"<total_steps>.*?</total_steps>",
            f"<total_steps> {total_steps} </total_steps>",
            text,
            flags=re.S,
        )
    return text


def run_ipi(src_xml, address, **overrides):
    """Run one i-PI + driver simulation and return the path to its ``data.out``.

    ``overrides`` may set ``temp`` (K), ``nbeads``, ``seed`` or ``total_steps``;
    anything left out keeps the template value. Each run happens in its own
    temporary directory so several can run concurrently without clashing.
    """
    tmp = tempfile.mkdtemp(prefix="ipi_run_")
    with open(os.path.join(tmp, "input.xml"), "w") as fh:
        fh.write(_prepare_input(src_xml, address, **overrides))

    sock = f"/tmp/ipi_{address}"
    if os.path.exists(sock):
        os.remove(sock)
    # Cap BLAS/OpenMP to one thread per child: we parallelise over trajectories
    # ourselves, so a multithreaded numpy/BLAS would oversubscribe the CPU.
    env = {
        **os.environ,
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }
    ipi = subprocess.Popen(
        ["i-pi", "input.xml"],
        cwd=tmp,
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    for _ in range(120):  # wait for the socket before launching the driver
        if os.path.exists(sock) or ipi.poll() is not None:
            break
        time.sleep(0.5)
    drv = subprocess.Popen(
        driver_command(address),
        cwd=tmp,
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    ipi.wait()
    drv.wait()
    return os.path.join(tmp, "data.out")


def temperature_for(bhw):
    """Kelvin temperature corresponding to a given beta*hbar*omega0."""
    return analysis.omega0() / (bhw * analysis.KELVIN_TO_HARTREE)


# %%
# Bosons and the energy-temperature curve
# ---------------------------------------
#
# We start with the headline feature: **bosons**. Listing all three atoms in
# ``<bosons> [0, 1, 2] </bosons>`` turns on bosonic exchange -- nothing else in
# the input changes. Exchange lets the three ring polymers connect into longer
# rings, which *lowers* the energy relative to distinguishable particles.
#
# Rather than a single number, let us trace the whole **energy-temperature
# curve** and compare it to the exact result. We sweep the dimensionless inverse
# temperature :math:`\beta\hbar\omega_0` from warm (near-classical) to cold,
# where the system settles into its quantum ground state.
#
# .. note::
#    **How many beads?** The number of beads needed grows with
#    :math:`\beta\hbar\omega_0`: the finite-bead (Trotter) error scales like
#    :math:`(\beta\hbar\omega_0/P)^2`, so a warm run needs far fewer beads than
#    a cold one. We therefore *scale* :math:`P` with the inverse temperature,
#    keeping every point converged to :math:`\sim0.2\%` -- and the warm point
#    runs with only 8 beads instead of 32.

# beta*hbar*omega0 and the (temperature-scaled) number of beads for each point
SWEEP_BHW = [1, 2, 3, 5]
SWEEP_BEADS = [8, 16, 24, 32]
SWEEP_STEPS = 6000  # short runs -- qualitative, reported with honest error bars
SKIP = 1000  # rows discarded as thermalisation (properties written every step)

omega0 = analysis.omega0()


def boson_point(bhw, nbeads):
    """Run bosons at inverse temperature ``bhw`` and return the mean total energy
    and its block-averaged statistical error, in units of hbar*omega0."""
    out = run_ipi(
        "input_3bosons.xml",
        f"bf-bsweep-{bhw}",
        temp=temperature_for(bhw),
        nbeads=nbeads,
        total_steps=SWEEP_STEPS,
    )
    mean, err = analysis.block_average(analysis.total_energy_series(out, SKIP))
    shutil.rmtree(os.path.dirname(out), ignore_errors=True)
    return mean / omega0, err / omega0


# The four temperatures are independent runs, so we launch them in parallel.
with concurrent.futures.ThreadPoolExecutor(max_workers=len(SWEEP_BHW)) as pool:
    sweep = list(
        pool.map(
            lambda i: boson_point(SWEEP_BHW[i], SWEEP_BEADS[i]),
            range(len(SWEEP_BHW)),
        )
    )
sweep_e = [m for m, _ in sweep]
sweep_err = [e for _, e in sweep]
for bhw, nbeads, (m, e) in zip(SWEEP_BHW, SWEEP_BEADS, sweep):
    print(
        f"bosons: beta*hbar*omega0 = {bhw}  (T = {temperature_for(bhw):5.1f} K, "
        f"P = {nbeads:2d})  ->  E/hw0 = {m:.3f} +/- {e:.3f}"
    )

# %%
# The PIMD points (with block-averaged error bars) sit on the exact bosonic curve
# (solid line). The distinguishable curve (dashed) is shown for reference:
# bosonic exchange always lies below it. At low temperature both approach the
# **ground-state energy** :math:`4.5\,\hbar\omega_0` (three particles, each
# contributing :math:`\tfrac{3}{2}\hbar\omega_0` of zero-point energy in 3D) --
# by :math:`\beta\hbar\omega_0 = 5` the boson energy is within ~0.5% of it. The
# runs are short, so the error bars are sizeable; they are what makes a scattered
# point "make sense" rather than look like a bug.

bhw_fine = np.linspace(0.5, 5.5, 80)
temps_fine = [temperature_for(b) for b in bhw_fine]

fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
ax.plot(
    bhw_fine,
    [analysis.analytical_energy(T, "bosonic") / omega0 for T in temps_fine],
    "-",
    color="C0",
    label="bosons (exact)",
)
ax.plot(
    bhw_fine,
    [analysis.analytical_energy(T, "dist") / omega0 for T in temps_fine],
    "--",
    color="gray",
    alpha=0.7,
    label="distinguishable (exact)",
)
ax.axhline(4.5, ls=":", color="k", alpha=0.5)
ax.text(
    0.55, 4.62, r"ground state $4.5\,\hbar\omega_0$", fontsize=8, color="k", alpha=0.7
)
ax.errorbar(
    SWEEP_BHW,
    sweep_e,
    yerr=sweep_err,
    fmt="o",
    color="C0",
    ms=8,
    capsize=4,
    label="bosons (PIMD)",
)
ax.set_xlabel(r"$\beta\hbar\omega_0$")
ax.set_ylabel(r"$\langle E\rangle\,/\,\hbar\omega_0$")
ax.set_title("Three bosons in a harmonic trap")
ax.legend()


# %%
# Switching statistics at a single temperature
# --------------------------------------------
#
# The whole point of the ``<bosons>`` tag is how *little* changes between cases.
# From the very same input we now compare, at one temperature
# (:math:`\beta\hbar\omega_0 = 2`, i.e. 17.4 K):
#
# * **distinguishable** particles -- ``<bosons> []`` (empty list),
# * **three bosons** -- ``<bosons> [0, 1, 2]``,
# * a **mixture** -- ``<bosons> [0, 1]``: atoms 0 and 1 exchange, atom 2 does not.
#
# Each is a separate template in ``data/`` differing only in that one line. We
# run all three at the same temperature and number of beads and read off the
# total energy.
#
# .. warning::
#    The **centroid-virial** kinetic estimator is *not* valid under bosonic
#    exchange, so the boson and mixture inputs record the primitive/quantum
#    virial (``virial_fq``) and thermodynamic (``kinetic_td``) estimators
#    instead. The total energy from the (always-valid) primitive virial is what
#    we compare below.

BHW_SWITCH = 2
BEADS_SWITCH = 16
T_SWITCH = temperature_for(BHW_SWITCH)

cases = [
    ("distinguishable", "input_3dist.xml", "dist", "bf-3dist"),
    ("bosons", "input_3bosons.xml", "bosonic", "bf-3bosons"),
    ("2 bosons + 1 dist", "input_2bosons1dist.xml", "mixed", "bf-2bosons1dist"),
]


def switch_point(xml, address):
    """Run one statistics case at the comparison temperature; return the mean
    total energy and its block-averaged error (Ha)."""
    out = run_ipi(
        xml, address, temp=T_SWITCH, nbeads=BEADS_SWITCH, total_steps=SWEEP_STEPS
    )
    mean, err = analysis.block_average(analysis.total_energy_series(out, SKIP))
    shutil.rmtree(os.path.dirname(out), ignore_errors=True)
    return mean, err


# the three cases are independent, so run them in parallel as well
with concurrent.futures.ThreadPoolExecutor(max_workers=len(cases)) as pool:
    switch = list(pool.map(lambda c: switch_point(c[1], c[3]), cases))
switch_sim = [m for m, _ in switch]
switch_err = [e for _, e in switch]
switch_ref = [analysis.analytical_energy(T_SWITCH, c[2]) for c in cases]
for (name, _, _, _), (m, e), ref in zip(cases, switch, switch_ref):
    print(
        f"{name:20s} (T = {T_SWITCH:.1f} K): PIMD {m * 1e3:.3f} +/- {e * 1e3:.3f} mHa  "
        f"exact {ref * 1e3:.3f} mHa"
    )

# %%
# Exchange orders the energies: at fixed temperature the bosonic energy is the
# lowest, distinguishable is higher, and the mixture sits in between -- each PIMD
# bar agreeing with its exact value within the (short-run) error bar.

labels = ["dist", "bosons", "2B+1D"]
x = np.arange(len(labels))
fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
ax.bar(
    x - 0.2,
    np.array(switch_sim) * 1e3,
    0.4,
    yerr=np.array(switch_err) * 1e3,
    capsize=4,
    label="PIMD",
)
ax.bar(x + 0.2, np.array(switch_ref) * 1e3, 0.4, label="exact")
ax.set_xticks(x, labels)
ax.set_ylabel("total energy / mHa")
ax.set_title(r"Three particles at $\beta\hbar\omega_0 = 2$ (17.4 K)")
ax.legend()


# %%
# Fermions: reweighting and the sign problem
# ------------------------------------------
#
# Fermions are harder. We cannot sample the antisymmetric density directly;
# instead we simulate **bosons** and reweight each sample by the fermionic sign
# :math:`s`,
#
# .. math::
#
#    \langle A \rangle_F = \frac{\langle A\, s \rangle}{\langle s \rangle},
#
# where the sign follows from the recursive configuration weight
#
# .. math::
#
#    W^{(N)} = \frac{1}{N}\sum_{k=1}^{N} \xi^{k-1}
#              e^{-\beta E_N^{(k)}} \, W^{(N-k)}, \qquad W^{(0)} = 1,
#
# with :math:`\xi = -1` for fermions. In the 2023 tutorial this weight was
# evaluated by a hand-written module; **i-PI 3.x computes it for us** and writes
# the sign to the ``fermionic_sign`` column, so we just request that property
# and reweight in two lines (see ``analysis.reweighted_fermionic_energy``).
#
# When the average sign :math:`\langle s\rangle` approaches zero the estimator
# becomes a ratio of two tiny, noisy numbers -- the fermionic **sign problem**.
# This is why the fermion input uses a higher temperature
# (:math:`\beta\hbar\omega_0 = 1.16`, 30 K) and fewer beads (12).
#
# .. admonition:: The exact reference value
#
#    The exact three-fermion energy at 30 K is **1.053 mHa**, computed by
#    ``analysis.analytical_energy`` from the canonical partition-function
#    recursion (elementary symmetric polynomial, :math:`\xi=-1`). An earlier
#    version of this tutorial used an incorrect hard-coded closed form that gave
#    0.912 mHa; the value here is confirmed independently by brute-force
#    enumeration of the three-fermion states. Note the Pauli exclusion principle
#    lifts the fermionic ground state to :math:`6.5\,\hbar\omega_0`, well above
#    the bosonic :math:`4.5\,\hbar\omega_0`.

fer_out = run_ipi("input_3fermions.xml", "bf-3fermions")
mean_sign, fer_energy = analysis.reweighted_fermionic_energy(fer_out, SKIP)
fer_ref = analysis.analytical_energy(30.0, "fermionic")
shutil.rmtree(os.path.dirname(fer_out), ignore_errors=True)

print("Three fermions (T = 30 K, single short run)")
print(f"  average sign <s>                  : {mean_sign:.4f}")
print(f"  fermionic total energy            : {fer_energy * 1e3:.4f} mHa")
print(f"  analytical total energy           : {fer_ref * 1e3:.4f} mHa")


# %%
# Don't trust one fermionic run -- average several
# ------------------------------------------------
#
# That single fermionic number is almost meaningless: the average sign is small,
# so the reweighted energy has a huge variance and one short run can land far
# from the truth. The honest thing is to run **several independent trajectories**
# (different random seeds) and combine them. Because the sign varies between
# trajectories, we use the sign-*weighted* estimator explained in the "Error
# estimation" section below (not a plain average). Here we run a handful of short
# trajectories in parallel -- still only a couple of minutes on a laptop.


def run_fermion_ensemble(n_traj, skip_rows):
    """Run ``n_traj`` short fermionic trajectories (different seeds) in parallel
    and combine them with the sign-weighted estimator. Returns
    ``(mean, error, n_eff, mean_sign, per_traj_energies)``."""
    seeds = [4001 + 137 * i for i in range(n_traj)]
    workers = min(n_traj, max(1, (os.cpu_count() or 2) // 2))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        outs = list(
            pool.map(
                lambda i: run_ipi(
                    "input_3fermions.xml", f"bf-fmulti-{i}", seed=seeds[i]
                ),
                range(n_traj),
            )
        )
    E, W, signs = [], [], []
    for out in outs:
        e_j, w_j, s_j = analysis.fermionic_trajectory_estimate(out, skip_rows)
        E.append(e_j)
        W.append(w_j)
        signs.append(s_j)
        shutil.rmtree(os.path.dirname(out), ignore_errors=True)
    mean, error, n_eff = analysis.weighted_average(E, W)
    return mean, error, n_eff, float(np.mean(signs)), np.array(E)


fer_mean, fer_err, n_eff, ens_sign, fer_traj = run_fermion_ensemble(
    n_traj=8, skip_rows=SKIP
)

print("Three fermions (T = 30 K, 8 short trajectories, sign-weighted)")
print(f"  average sign <s>                  : {ens_sign:.3f}")
print(
    f"  weighted fermionic energy         : {fer_mean * 1e3:.4f} +/- {fer_err * 1e3:.4f} mHa"
)
print(f"  analytical total energy           : {fer_ref * 1e3:.4f} mHa")
print(f"  (effective sample size n_eff = {n_eff:.1f} of 8)")

# %%
# The individual trajectories (grey points) scatter wildly -- some land far from
# the exact value -- yet the sign-weighted average (blue, with its error bar)
# *brackets* the analytical result. That large-but-honest error bar is the
# fingerprint of the fermionic **sign problem**: the numbers make sense, they are
# just noisy, and tightening them needs far more sampling. Contrast this with the
# distinguishable/bosonic runs above, whose single short runs were already close
# because there is no sign to fight.

fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
ax.axhline(fer_ref * 1e3, color="k", ls="--", label="exact (1.053 mHa)")
ax.plot(
    np.zeros_like(fer_traj) + 0.15 * (np.arange(len(fer_traj)) - 3.5) / 3.5,
    fer_traj * 1e3,
    "o",
    color="gray",
    alpha=0.6,
    label="individual trajectories",
)
ax.errorbar(
    [0],
    [fer_mean * 1e3],
    yerr=[fer_err * 1e3],
    fmt="s",
    color="C0",
    ms=9,
    capsize=5,
    label="sign-weighted mean",
)
ax.set_xlim(-1, 1)
ax.set_xticks([])
ax.set_ylabel("fermionic total energy / mHa")
ax.set_title("Three fermions at 30 K: one run vs. eight")
ax.legend()


# %%
# Error estimation for fermions: weight trajectories by the sign
# --------------------------------------------------------------
#
# For fermions a *single* run gives one number with no error bar, and naively
# averaging the per-trajectory reweighted energies
# :math:`E_j = \langle \varepsilon\, s\rangle_j / \langle s\rangle_j` with a
# plain :math:`\mathrm{std}/\sqrt{M}` is **biased**: a trajectory that happened
# to sample mostly low-weight (small-sign) configurations yields a wild ratio,
# yet counts as much as a well-converged one.
#
# The fix (SI of Hirshberg, Invernizzi & Parrinello, *J. Chem. Phys.* **152**,
# 171102 (2020)) is to weight each trajectory by its total weight in the
# reweighted ensemble, :math:`W_j = \sum_\mathrm{steps} s` (the sum of the
# instantaneous signs), and to use an *effective* sample size:
#
# .. math::
#
#    \bar E_F = \frac{\sum_j W_j E_j}{\sum_j W_j}, \qquad
#    n_\mathrm{eff} = \frac{\left(\sum_j W_j\right)^2}{\sum_j W_j^2}, \qquad
#    \sigma_E^2 = \frac{n_\mathrm{eff}}{n_\mathrm{eff}-1}
#                 \frac{\sum_j W_j (E_j - \bar E_F)^2}{\sum_j W_j},
#
# with statistical error :math:`\sigma_E/\sqrt{n_\mathrm{eff}}`. When all
# weights are equal this reduces to the ordinary mean and
# :math:`\mathrm{std}/\sqrt{M}`, which is why the bosonic runs above needed no
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
# naive mean-of-ratios        1.45  (biased high)    :math:`\pm 0.32`
# SI weighted, assuming n=M    1.11                   :math:`\pm 0.078`
# SI weighted, using n_eff     1.11                   :math:`\pm 0.087`  (n_eff=16 of 20)
# ==========================  =====================  =============================
#
# Two lessons: the weighting **removes the bias** in the mean (the naive
# mean-of-ratios sits at 1.45, the weighted estimate at 1.11), and using the
# effective sample size gives an **honest, ~12% larger** error bar
# (:math:`\sqrt{M/n_\mathrm{eff}}`) than pretending all :math:`M` trajectories
# are equally informative. The weighted estimate 1.11 :math:`\pm` 0.09 agrees
# with the exact 1.053 mHa within its error bar -- there is *no* mysterious
# fermionic discrepancy once (i) the correct benchmark is used and (ii) the
# statistical error is estimated properly.
#
# The runs in this tutorial are deliberately short and use only 12 beads for the
# fermions -- enough to "make sense", not to be tightly converged. For production
# you would simply run more (and longer) trajectories and, if needed, increase
# the number of beads; the estimators above are exactly what you would use.


# %%
#
# **References**
#
# * B. Hirshberg, V. Rizzi, M. Parrinello,
#   `PNAS 116, 21445 (2019) <https://doi.org/10.1073/pnas.1913365116>`_
# * Y. M. Y. Feldman, B. Hirshberg,
#   `JCP 159, 154107 (2023) <https://doi.org/10.1063/5.0173749>`_
# * B. Hirshberg, M. Invernizzi, M. Parrinello,
#   `JCP 152, 171102 (2020) <https://doi.org/10.1063/5.0008720>`_
