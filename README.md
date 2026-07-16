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
  potential, so there is no Fortran to build. (If the compiled `i-pi-driver` is
  on your `PATH`, the recipe auto-detects and uses it — a few times faster.)
* **Native fermions.** i-PI 3.x records the `fermionic_sign` directly, so the
  hand-written reweighting module (and its `MDAnalysis` dependency) is gone —
  fermionic averages are two lines of NumPy.

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
problem** — where a single short run is nearly useless and you must average
several trajectories with a proper (sign-weighted) error bar.

## Quick start (conda — recommended for the tutorial)

```bash
git clone https://github.com/BarakHirshberg/PIQM2026_BH_Tutorial && cd PIQM2026_BH_Tutorial

# create a local environment (a folder ./env, keeps your base conda clean)
conda env create --prefix ./env --file examples/bosons-fermions-pimd/environment.yml
conda activate ./env

# run the whole recipe (a few minutes on a laptop)
cd examples/bosons-fermions-pimd
python bosons-fermions-pimd.py
```

The script runs the four cases, prints each PIMD energy next to the **exact**
value, runs a small fermionic multi-trajectory average, and saves the figures.

### Run it as a notebook

The recipe is in the sphinx-gallery "percent" format: every `# %%` marks a cell,
so it opens directly as a notebook in **VS Code** or **Jupyter** (via
`jupytext`):

```bash
conda activate ./env
pip install jupytext jupyterlab
cd examples/bosons-fermions-pimd
jupytext --to notebook bosons-fermions-pimd.py
jupyter lab bosons-fermions-pimd.ipynb
```

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
├── analysis.py                 # output reader + EXACT energies + fermionic reweighting + weighted error
└── data/                       # the four i-PI input files
```

The example directory follows the
[atomistic-cookbook](https://github.com/lab-cosmo/atomistic-cookbook)
conventions (a `README.rst`, an `environment.yml`, and a sphinx-gallery `.py`
that renders to a notebook + HTML page), so it can be contributed upstream by
copying it into `examples/` of that repository.

## The exact benchmarks (and a fixed bug)

Each case is compared against the **exact** energy of non-interacting particles
in a harmonic trap, computed in `analysis.py` from the canonical
partition-function recursion (mHa):

| Case | exact energy |
|------|--------------|
| 3 distinguishable (17.4 K) | 0.6514 |
| 3 bosons (17.4 K) | 0.5803 |
| 2 bosons + 1 dist (17.4 K) | 0.6235 |
| 3 fermions (30 K) | **1.0530** |

> **Note.** The 2023 tutorial hard-coded an *incorrect* closed form for the
> three-fermion energy (0.912 mHa). The correct value is **1.053 mHa**, verified
> two independent ways (canonical recursion + brute-force enumeration of the
> three-fermion states). This wrong benchmark — not the simulation — was the
> main reason fermions previously appeared to "disagree." With the correct
> value the PIMD result agrees within its error bar.

## Fermions need averaging and careful error bars

A single short fermionic run is nearly meaningless: the average sign is small
(⟨s⟩ ≈ 0.3 — the **sign problem**), so the reweighted energy has a huge
variance. The tutorial runs a handful of short trajectories and combines them
with the **sign-weighted** estimator from the SI of Hirshberg, Invernizzi &
Parrinello (*JCP* **152**, 171102, 2020): weight each trajectory by its total
sign `W_j = Σ s`, use an effective sample size `n_eff = (Σ W)² / Σ W²`, and

```
Ē_F = Σ Wⱼ Eⱼ / Σ Wⱼ ,   σ² = n/(n−1) · Σ Wⱼ(Eⱼ−Ē_F)²/Σ Wⱼ ,   error = σ/√n_eff
```

Naively averaging the per-trajectory ratios with `std/√M` is *biased* (a
low-weight trajectory gives a wild ratio yet counts equally). Post-processing 20
trajectories × 5000 steps (exact 1.053 mHa):

| estimator | mean (mHa) | error (mHa) |
|-----------|-----------|-------------|
| naive mean-of-ratios | 1.45 (biased high) | ± 0.32 |
| SI weighted (using n_eff) | 1.11 | ± 0.09 (n_eff = 16 of 20) |

The weighting removes the bias, and the weighted estimate **1.11 ± 0.09 agrees
with the exact 1.053** within error. `weighted_average()` in `analysis.py`
reduces to the plain mean + `std/√M` when weights are equal, so the same code
handles bosons and the average sign.

The tutorial runs are deliberately short and use only 12 beads for fermions —
enough to *make sense*, not to be tightly converged. For production you would
simply run more (and longer) trajectories, and increase the number of beads if
needed; the estimators are the same. Building the compiled **f90 driver**
(`i-pi-driver -m harm3d`, `make` in the i-PI source `drivers/f90`) makes such
runs several times faster, and the recipe uses it automatically if it is on
your `PATH`.

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
