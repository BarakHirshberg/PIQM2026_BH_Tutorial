"""
Path integral molecular dynamics of bosons and fermions
========================================================

:Authors: Barak Hirshberg `@bhirshberg <https://github.com/bhirshberg>`_

This example shows how to run path integral molecular dynamics (PIMD)
simulations of **indistinguishable** particles -- bosons and fermions -- with
``i-PI``, and how to analyze the output. It is a modernised version of the
"Bosonic and Fermionic PIMD" module of the `PIQM 2023 tutorial
<https://github.com/i-pi/piqm2023-tutorial>`_, updated to i-PI 3.x.

We simulate a few non-interacting particles (mass 1) in a three-dimensional
isotropic harmonic trap with :math:`\\hbar\\omega_0 = 3\\,\\mathrm{meV}` (a very
soft trap) and study how quantum statistics changes the average energy. Most of
the tutorial uses three particles; the statistics-comparison section uses four,
so the energy differences between the cases are large enough to see clearly. We
work in the natural dimensionless inverse temperature
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

import numpy as np

import analysis
from scripts import plots


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
# How we run each case
# --------------------
#
# The mechanics of driving i-PI from Python -- launching the ``i-pi`` server and
# a force driver that talk over a socket, and editing the XML input -- are not
# the subject of this tutorial, so they live in ``scripts/ipi_runs.py``. From
# here we only need two calls:
#
# .. code-block:: python
#
#     out = run_ipi("input_3bosons.xml", temp=17.4, nbeads=16)   # -> path to data.out
#     outs = run_parallel(jobs)                                  # several at once
#
# ``run_ipi`` starts from one of the template inputs in ``data/`` (one per case),
# overrides only the knobs we vary -- which atoms exchange, temperature, beads,
# seed, length -- runs the simulation with i-PI's built-in ``harmonic`` driver,
# and returns the output file for analysis. Every case in this tutorial runs
# several inputs at once, so we import ``run_parallel``; everything below is
# physics and analysis.

from scripts.ipi_runs import run_parallel


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

# Run bosons at each temperature: the *same* input (``<bosons> [0, 1, 2]``), only
# the temperature and bead count change. The runs are independent, so we launch
# them together.
jobs = [
    (
        "input_3bosons.xml",
        f"bsweep-{bhw}",
        dict(temp=analysis.temperature_for(bhw), nbeads=P, total_steps=SWEEP_STEPS),
    )
    for bhw, P in zip(SWEEP_BHW, SWEEP_BEADS)
]
outputs = run_parallel(jobs)

# Analysis: block-averaged total energy (the estimator valid under exchange),
# in units of hbar*omega0.
sweep = [analysis.block_average(analysis.total_energy_series(o, SKIP)) for o in outputs]
sweep_e = [m / omega0 for m, _ in sweep]
sweep_err = [e / omega0 for _, e in sweep]
for bhw, P, m, e in zip(SWEEP_BHW, SWEEP_BEADS, sweep_e, sweep_err):
    print(
        f"beta*hbar*omega0 = {bhw}  (T = {analysis.temperature_for(bhw):5.1f} K, "
        f"P = {P:2d})  ->  E/hw0 = {m:.3f} +/- {e:.3f}"
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

plots.plot_boson_energy_curve(SWEEP_BHW, sweep_e, sweep_err)


# %%
# Switching statistics at a single temperature
# --------------------------------------------
#
# The whole point of the ``<bosons>`` tag is how *little* changes between cases.
# We now compare, at one temperature (:math:`\beta\hbar\omega_0 = 2`, i.e.
# 17.4 K), **four** particles under three kinds of statistics. We use four here
# (rather than the three of the other sections) so the exchange effect is large
# enough that the energies separate clearly given these short runs:
#
# * **distinguishable** -- ``<bosons> []`` (empty list),
# * **four bosons** -- ``<bosons> [0, 1, 2, 3]``,
# * a **mixture** -- ``<bosons> [0, 1, 2]``: atoms 0-2 exchange, atom 3 does not.
#
# Each is a separate template in ``data/`` differing only in that one line:
#
# .. code-block:: xml
#
#     <bosons> []           </bosons>   <!-- distinguishable: no exchange -->
#     <bosons> [0, 1, 2, 3] </bosons>   <!-- four bosons: every pair may exchange -->
#     <bosons> [0, 1, 2]    </bosons>   <!-- mixture: 0-2 exchange, atom 3 stays distinct -->
#
# Exchange means the ring polymers of the listed atoms are allowed to connect
# into longer rings; the atoms left out stay closed on themselves.
#
# .. warning::
#    The **centroid-virial** kinetic estimator is *not* valid under bosonic
#    exchange, so the boson and mixture inputs record the primitive/quantum
#    virial (``virial_fq``) and thermodynamic (``kinetic_td``) estimators
#    instead. The total energy from the (always-valid) primitive virial is what
#    we compare below.
#
# You can see this in the ``<properties>`` line of each input: the
# distinguishable case keeps ``kinetic_cv``, while the exchange cases drop it.
#
# .. code-block:: xml
#
#     <!-- distinguishable -->
#     <properties ...> [ ..., kinetic_cv, kinetic_td, potential, virial_fq ] </properties>
#
#     <!-- bosons / mixture: no kinetic_cv -->
#     <properties ...> [ ..., kinetic_td, potential, virial_fq ] </properties>

BHW_SWITCH = 2
BEADS_SWITCH = 16
T_SWITCH = analysis.temperature_for(BHW_SWITCH)

# (label, input template, (n_bosons, n_dist) for the exact reference)
cases = [
    ("distinguishable", "input_4dist.xml", (0, 4)),
    ("4 bosons", "input_4bosons.xml", (4, 0)),
    ("3 bosons + 1 dist", "input_3bosons1dist.xml", (3, 1)),
]
jobs = [
    (
        xml,
        f"switch-{i}",
        dict(temp=T_SWITCH, nbeads=BEADS_SWITCH, total_steps=SWEEP_STEPS),
    )
    for i, (_, xml, _) in enumerate(cases)
]
outputs = run_parallel(jobs)

switch = [
    analysis.block_average(analysis.total_energy_series(o, SKIP)) for o in outputs
]
switch_sim = [m for m, _ in switch]
switch_err = [e for _, e in switch]
switch_ref = [analysis.mixture_energy(nb, nd, T_SWITCH) for _, _, (nb, nd) in cases]
for (name, _, _), m, e, ref in zip(cases, switch_sim, switch_err, switch_ref):
    print(
        f"{name:20s} (T = {T_SWITCH:.1f} K): PIMD {m * 1e3:.3f} +/- {e * 1e3:.3f} mHa  "
        f"exact {ref * 1e3:.3f} mHa"
    )

# %%
# Exchange orders the energies: at fixed temperature the bosonic energy is the
# lowest, distinguishable is highest, and the mixture sits in between. With four
# particles the gaps are large enough to resolve with these short runs, and each
# PIMD bar matches its exact value.

plots.plot_statistics_bars(
    ["dist", "4 bosons", "3B+1D"],
    switch_sim,
    switch_err,
    switch_ref,
    r"Four particles at $\beta\hbar\omega_0 = 2$ (17.4 K)",
)


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
# In the input, there is no separate "fermion mode": you run an ordinary
# *bosonic* simulation and simply add ``fermionic_sign`` to the property list:
#
# .. code-block:: xml
#
#     <bosons> [0, 1, 2] </bosons>                     <!-- still a bosonic run -->
#     <properties ...> [ ..., virial_fq, fermionic_sign ] </properties>
#
# When the average sign :math:`\langle s\rangle` approaches zero the estimator
# becomes a ratio of two tiny, noisy numbers -- the fermionic **sign problem**.
# This is why the fermion input uses a higher temperature
# (:math:`\beta\hbar\omega_0 = 1.16`, 30 K) and fewer beads (12).
#
# .. note::
#    **The exact reference value.** The exact three-fermion energy at 30 K is
#    **1.053 mHa**, computed by ``analysis.analytical_energy`` from the canonical
#    partition-function recursion (elementary symmetric polynomial,
#    :math:`\xi=-1`). An earlier version of this tutorial used an incorrect
#    hard-coded closed form that gave 0.912 mHa; the value here is confirmed
#    independently by brute-force enumeration of the three-fermion states. Note
#    the Pauli exclusion principle lifts the fermionic ground state to
#    :math:`6.5\,\hbar\omega_0`, well above the bosonic :math:`4.5\,\hbar\omega_0`.
#
# We estimate the fermionic energy from **eight independent trajectories**
# (different random seeds), combined into a single sign-weighted average with an
# error bar (the reason for the sign-weighting is spelled out afterwards).

# Eight independent fermionic trajectories with different random seeds, together.
seeds = [4001 + 137 * i for i in range(8)]
jobs = [
    ("input_3fermions.xml", f"fmulti-{i}", dict(seed=s)) for i, s in enumerate(seeds)
]
outputs = run_parallel(jobs)

# Per-trajectory reweighted energy E_j and total sign weight W_j, combined with
# the sign-weighted estimator (see "How the error bar is estimated" below).
E, W, signs = [], [], []
for o in outputs:
    e_j, w_j, s_j = analysis.fermionic_trajectory_estimate(o, SKIP)
    E.append(e_j)
    W.append(w_j)
    signs.append(s_j)
fer_mean, fer_err, n_eff = analysis.weighted_average(E, W)
fer_traj = np.array(E)
fer_ref = analysis.analytical_energy(30.0, "fermionic")

print("Three fermions (T = 30 K, 8 short trajectories)")
print(f"  average sign <s>          : {np.mean(signs):.3f}")
print(f"  fermionic energy          : {fer_mean * 1e3:.3f} +/- {fer_err * 1e3:.3f} mHa")
print(f"  exact                     : {fer_ref * 1e3:.3f} mHa")
print(f"  effective sample size     : n_eff = {n_eff:.1f} of 8")

# %%
# The 8-trajectory fermionic energy brackets the exact 1.053 mHa within its error
# bar. The individual trajectories (grey) scatter around the sign-weighted mean
# (blue, with its error bar), which sits on the exact value -- the large-but-honest
# error bar is the fingerprint of the fermionic **sign problem**.

plots.plot_fermion_ensemble(fer_traj, fer_mean, fer_err, fer_ref)


# %%
# How the error bar is estimated
# ------------------------------
#
# Combining fermionic trajectories takes more care than a plain average, because
# they carry different *effective* amounts of information: a trajectory that
# sampled a smaller average sign has a more poorly-determined ratio and must count
# for less. We therefore weight each trajectory :math:`j` by its total sign
# :math:`W_j = \sum_\mathrm{steps} s`, and measure an *effective sample size*:
#
# .. math::
#
#    \bar E_F = \frac{\sum_j W_j E_j}{\sum_j W_j}, \qquad
#    n_\mathrm{eff} = \frac{\left(\sum_j W_j\right)^2}{\sum_j W_j^2}, \qquad
#    \sigma_E^2 = \frac{n_\mathrm{eff}}{n_\mathrm{eff}-1}
#                 \frac{\sum_j W_j (E_j - \bar E_F)^2}{\sum_j W_j},
#
# with statistical error :math:`\sigma_E/\sqrt{n_\mathrm{eff}}`
# (SI of Hirshberg, Invernizzi & Parrinello, *J. Chem. Phys.* 152, 171102, 2020).
# The weighted mean is exactly what you would get by pooling every sample of every
# trajectory into one long run, :math:`\sum \varepsilon s / \sum s`. The effective
# sample size came out :math:`n_\mathrm{eff} \approx 7` of 8 above: the
# trajectories are worth fewer than eight fully-independent samples because they
# carry unequal weight, so the honest error bar is a little larger than
# :math:`\mathrm{std}/\sqrt{8}` would give. When all the weights are equal -- as
# for the sign-free bosons -- :math:`n_\mathrm{eff} = M` and this reduces to the
# ordinary mean and standard error, which is why the boson runs needed no special
# treatment. Both estimators are in ``analysis.py`` (:func:`weighted_average`,
# :func:`fermionic_trajectory_estimate`).


# %%
# The sign problem, made concrete
# -------------------------------
#
# The sign problem is not abstract: lower the temperature and the average sign
# shrinks, so the *same* eight trajectories scatter far more widely. Below we
# repeat the run at 20 K and put the per-trajectory energies next to the 30 K ones.

cold_jobs = [
    ("input_3fermions.xml", f"fcold-{i}", dict(temp=20.0, seed=s))
    for i, s in enumerate(seeds)
]
cold_E, cold_signs = [], []
for o in run_parallel(cold_jobs):
    e_j, _, s_j = analysis.fermionic_trajectory_estimate(o, SKIP)
    cold_E.append(e_j)
    cold_signs.append(s_j)
cold_traj = np.array(cold_E)
cold_ref = analysis.analytical_energy(20.0, "fermionic")

print(
    f"average sign:  30 K -> {np.mean(signs):.3f},   20 K -> {np.mean(cold_signs):.3f}"
)

plots.plot_sign_scatter(
    [
        (f"30 K\n<s> = {np.mean(signs):.2f}", fer_traj, fer_ref),
        (f"20 K\n<s> = {np.mean(cold_signs):.2f}", cold_traj, cold_ref),
    ]
)

# %%
# At 20 K the average sign has dropped noticeably (from ~0.36 to ~0.24) and the
# per-trajectory energies scatter over a far wider range -- some even come out
# negative or well above the true value. Eight short trajectories no longer pin
# the energy down; you would
# need far more sampling. This is the sign problem, and it is why we run the
# fermions at the milder 30 K -- colder still, or with more particles, both the
# scatter and the number of trajectories you need grow exponentially.


# %%
# The runs here are deliberately short and use only 12 beads for the fermions --
# enough to "make sense", not to be tightly converged. For production you would
# run more (and longer) trajectories and, if needed, more beads; the sign-weighted
# estimator above is exactly what you would use.


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
