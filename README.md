# PIMD of bosons and fermions with i-PI — PIQM 2026 tutorial

Path integral molecular dynamics (PIMD) simulations of **indistinguishable
particles** (bosons and fermions) with [i-PI](https://ipi-code.org) 3.x,
prepared for the **PIQM 2026** meeting in Shanghai.

This is a modernised, self-contained rewrite of the "Bosonic and Fermionic
PIMD" module of the [PIQM 2023 tutorial](https://github.com/i-pi/piqm2023-tutorial),
laid out as an [Atomistic Cookbook](https://atomistic-cookbook.org) recipe.
Two things make it far easier to run than the 2023 version:

* **No compilation, no special branch.** Everything installs from PyPI. i-PI 3.x
  ships a pure-Python driver (`i-pi-py_driver`) with a built-in `harmonic`
  potential, so there is no Fortran to build.
* **Native fermions.** i-PI 3.x records the `fermionic_sign` directly, so the
  hand-written reweighting module (and its `MDAnalysis` dependency) is gone —
  fermionic averages are now two lines of NumPy.

## What you will do

Simulate three non-interacting particles in a 3D harmonic trap and compute the
average energy for four kinds of statistics, each a **one-line change** to the
`<bosons>` tag in the i-PI input:

| Case | `<bosons>` | beads | T (K) |
|------|-----------|-------|-------|
| 3 distinguishable | `[]` | 32 | 17.4 |
| 3 bosons | `[0, 1, 2]` | 32 | 17.4 |
| 2 bosons + 1 distinguishable | `[0, 1]` | 32 | 17.4 |
| 3 fermions (reweighted) | `[0, 1, 2]` + `fermionic_sign` | 12 | 30.0 |

You will see that exchange orders the energies as
E(bosons) < E(distinguishable) < E(fermions), and meet the fermionic **sign
problem**.

## Quick start (conda — recommended for the tutorial)

```bash
git clone <this-repo> && cd PIQM2026_BH_Tutorial

# create a local environment (a folder ./env, keeps your base conda clean)
conda env create --prefix ./env --file examples/bosons-fermions-pimd/environment.yml
conda activate ./env

# run the whole recipe (~2–3 minutes)
cd examples/bosons-fermions-pimd
python bosons-fermions-pimd.py
```

The script runs all four simulations and prints, for each, the PIMD energy next
to the exact analytical value, then saves two figures.

### Run it as a notebook

The recipe is written in the sphinx-gallery "percent" format: every `# %%`
marks a new cell, so the file opens directly as a notebook in **VS Code** or
**Jupyter** (via `jupytext`):

```bash
conda activate ./env
pip install jupytext jupyterlab
cd examples/bosons-fermions-pimd
jupytext --to notebook bosons-fermions-pimd.py
jupyter lab bosons-fermions-pimd.ipynb
```

(The narrative text is carried along as comments; it is rendered as formatted
markdown in the HTML page / downloadable notebook produced by the cookbook
build.)

### Pip-only alternative (no conda)

All dependencies are on PyPI, so a plain virtual environment works too:

```bash
python -m venv env && source env/bin/activate
pip install "ipi>=3.2.0" numpy matplotlib ase "chemiscope>=1.0"
cd examples/bosons-fermions-pimd && python bosons-fermions-pimd.py
```

Or, with [nox](https://nox.thea.codes) installed, `nox -e bosons-fermions-pimd`
builds the environment and runs the recipe for you.

## Repository layout

```
examples/bosons-fermions-pimd/
├── README.rst                  # one-line recipe description
├── environment.yml             # conda/pip dependencies
├── bosons-fermions-pimd.py     # the tutorial (sphinx-gallery format)
├── analysis.py                 # output reader + analytical references + fermionic reweighting
├── data/                       # the four i-PI input files
└── reference/
    └── run_convergence.py      # offline: N trajectories/case → mean ± error table
```

The example directory follows the
[atomistic-cookbook](https://github.com/lab-cosmo/atomistic-cookbook)
conventions (a `README.rst`, an `environment.yml`, and a sphinx-gallery `.py`
that renders to a notebook + HTML page), so it can be contributed upstream by
copying it into `examples/` of that repository.

## Convergence and error bars

A single short run is deliberately noisy. Because the dynamics are unbiased
(NVT with a thermostat), averaging several **independent** trajectories
converges toward the exact value and gives a proper standard error. The helper
`reference/run_convergence.py` runs this study in parallel across your CPU
cores:

```bash
cd examples/bosons-fermions-pimd
python reference/run_convergence.py          # 10 trajectories per case
```

With 10 trajectories (same short settings as the tutorial) the energies land
within about one standard error of the analytical values (mHa):

| Case | trajectory average | analytical |
|------|--------------------|------------|
| 3 distinguishable | 0.643 ± 0.017 | 0.651 |
| 3 bosons | 0.569 ± 0.017 | 0.580 |
| 2 bosons + 1 dist | 0.606 ± 0.018 | 0.624 |
| 3 fermions | 1.11 ± 0.09 (weighted) | 0.912 |

(10 trajectories per case.) The three distinguishable/bosonic cases land within
~1σ of the exact values — for equal-weight data, averaging over trajectories
just removes the noise of a single run.

### Error estimation for fermions

The **fermionic** estimator is fundamentally harder (small average sign ⟨s⟩ ≈
0.3 — the **sign problem**), and its error must be estimated carefully. Naively
averaging the per-trajectory reweighted energies with `std/√M` is *biased*: a
trajectory that sampled mostly low-weight configurations gives a wild ratio but
counts equally. Following the SI of Hirshberg, Invernizzi & Parrinello (*JCP*
**152**, 171102, 2020), each trajectory is weighted by its total sign
`W_j = Σ s`, and an effective sample size `n_eff = (Σ W)² / Σ W²` is used:

```
Ē_F = Σ Wⱼ Eⱼ / Σ Wⱼ ,   σ² = n/(n−1) · Σ Wⱼ(Eⱼ−Ē_F)²/Σ Wⱼ ,   error = σ/√n_eff
```

Post-processing 20 trajectories × 5000 steps (analytical 0.912 mHa):

| estimator | mean (mHa) | error (mHa) |
|-----------|-----------|-------------|
| naive mean-of-ratios | 1.45 (biased high) | ± 0.32 |
| SI weighted, assuming n = M | 1.11 | ± 0.078 |
| SI weighted, using n_eff | 1.11 | ± 0.087 (n_eff = 16 of 20) |

The weighting **removes the bias** in the mean, and using `n_eff` gives an
**honest, ~12% larger** error bar than assuming all M trajectories count fully.
`weighted_average()` in `analysis.py` reduces to the plain mean + `std/√M` when
weights are equal, so the same code handles bosons and the average sign. Run:

```bash
N_TRAJ=20 STEPS=5000 MAX_CONCURRENT=12 CASE_FILTER=fermion python reference/run_convergence.py
```

For production-scale fermionic runs the compiled **f90 driver** (`i-pi-driver
-m harm3d`, built from the i-PI source with `make` in `drivers/f90`) is several
times faster than the Python driver and is the practical choice; the tutorial
uses the Python driver only to stay pip-installable.

## Background and references

* Bosonic PIMD: B. Hirshberg, V. Rizzi, M. Parrinello,
  [*PNAS* **116**, 21445 (2019)](https://doi.org/10.1073/pnas.1913365116)
* Quadratic-scaling exchange: Y. M. Y. Feldman, B. Hirshberg,
  [*J. Chem. Phys.* **159**, 154107 (2023)](https://doi.org/10.1063/5.0173749)
* Fermions by reweighting: Hirshberg, Invernizzi, Parrinello,
  [*J. Chem. Phys.* **152**, 171102 (2020)](https://doi.org/10.1063/5.0008720)
* i-PI: [ipi-code.org](https://ipi-code.org)
* Original tutorial: [i-pi/piqm2023-tutorial](https://github.com/i-pi/piqm2023-tutorial)

## License

BSD 3-Clause (see [LICENSE](LICENSE)).
