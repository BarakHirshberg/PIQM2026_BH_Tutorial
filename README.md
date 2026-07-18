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
  on your `PATH`, the recipe auto-detects and uses it, though for this tiny
  three-atom system it makes little difference — i-PI's per-step overhead
  dominates, not the force evaluation.)
* **Native fermions.** i-PI 3.x records the `fermionic_sign` directly, so the
  hand-written reweighting module (and its `MDAnalysis` dependency) is gone —
  fermionic averages are two lines of NumPy.

## What you will do

Simulate a few non-interacting particles ($m=1$) — **three** in most sections,
**four** in the statistics comparison — in a **3D isotropic harmonic trap** with
force constant `k = 1.21647924e-8` Ha/Bohr², i.e. a trap frequency `ℏω₀ = 3 meV`
(a very soft trap). Quantum statistics is switched on with a **one-line change**
to the `<bosons>` tag in the i-PI input, and temperature is quoted in the natural
dimensionless unit `βℏω₀ = ℏω₀ / k_B T`. The tutorial is in three steps:

1. **Bosons and the energy–temperature curve.** Turn on exchange with
   `<bosons> [0, 1, 2]` and trace `⟨E⟩` across a sweep of temperatures,
   `βℏω₀ = 1, 2, 3, 5`, comparing the PIMD points against the exact curve and
   watching the energy settle onto the ground state at low `T`.
2. **Switching statistics.** From the *same* input flip to `[]`
   (distinguishable), `[0, 1, 2, 3]` (four bosons), or `[0, 1, 2]` (three bosons
   + one distinguishable) and compare the energies at one temperature
   (`βℏω₀ = 2`, i.e. 17.4 K). This section uses **four** particles so the exchange
   differences are large enough to resolve with short runs.
3. **Fermions.** Reweight the bosonic run by the `fermionic_sign`, meet the
   **sign problem**, and average several trajectories to get an honest
   (sign-weighted) error bar. Run warmer (`βℏω₀ = 1.16`, 30 K) with fewer beads.

Colder runs need more beads `P`, so instead of a fixed value we scale `P` with
`βℏω₀`. Here are the conditions (`βℏω₀`, temperature, beads) used in each part:

| `βℏω₀` | T (K) | beads `P` | role |
|--------|-------|-----------|------|
| 1 | 34.8 | 8 | boson sweep (warm) |
| 2 | 17.4 | 16 | boson sweep + statistics comparison |
| 3 | 11.6 | 24 | boson sweep |
| 5 | 7.0 | 32 | boson sweep (cold — near ground state) |
| 1.16 | 30.0 | 12 | fermions (sign-problem limited) |

You will see that exchange orders the energies as
E(bosons) < E(distinguishable) < E(fermions), and meet the fermionic **sign
problem** — which is why the fermions are estimated from several trajectories
combined with a proper (sign-weighted) error bar.

## Before the tutorial: set up and verify (please do this ahead of time)

Set up **before the meeting** — installing at the venue over shared WiFi is slow
and failure-prone. You need `conda` (Miniconda or Anaconda).

```bash
# 1. clone and enter the example folder
git clone https://github.com/BarakHirshberg/PIQM2026_BH_Tutorial
cd PIQM2026_BH_Tutorial/examples/bosons-fermions-pimd

# 2. create a local environment (a ./env folder; keeps your base conda clean)
conda env create --prefix ./env --file environment.yml
conda activate ./env

# 3. tools to read/run it as a notebook
pip install sphinx-gallery jupyterlab
```

**Verify (≈15 s).** This runs a tiny 50-step simulation end to end — the i-PI
server, the force driver, and their socket — and should print `setup OK`:

```bash
python -c "import analysis; from scripts.ipi_runs import run_ipi; \
o = run_ipi('input_3bosons.xml', total_steps=50, nbeads=8); \
print('setup OK' if analysis.read_ipi_output(o).get('potential') is not None else 'FAIL')"
```

If you see `setup OK`, you are ready. If it fails, email
<barak.hirshberg@gmail.com> with the error message **before** the meeting.

## During the tutorial: open the notebook

In a terminal, go to the example folder inside your clone (adjust the path if you
cloned somewhere other than your home directory) and activate the environment,
then convert the recipe to a Jupyter notebook (proper markdown + rendered math,
code in code cells) with sphinx-gallery's own converter and run the cells top to
bottom — the four figures appear inline:

```bash
cd PIQM2026_BH_Tutorial/examples/bosons-fermions-pimd   # the folder with your ./env
conda activate ./env
sphinx_gallery_py2jupyter bosons-fermions-pimd.py   # -> bosons-fermions-pimd.ipynb
jupyter lab bosons-fermions-pimd.ipynb              # then Run -> Run All Cells
```

The full run takes **~5–15 minutes** (the fermion trajectories are the slow
part). You can read all the text immediately; the figures fill in as the cells
finish.

**Headless or remote (no browser)?** Run the notebook non-interactively, or just
run the script for the printed numbers (figures are captured but not shown):

```bash
jupyter nbconvert --to notebook --execute bosons-fermions-pimd.ipynb   # runs every cell
python bosons-fermions-pimd.py                                         # numbers only
```

> Use the sphinx-gallery converter, **not** `jupytext --to notebook`: jupytext
> leaves the reStructuredText prose as raw comments inside code cells instead of
> rendered markdown. In **VS Code** you can instead open the `.py` directly (the
> `# %%` markers give "Run Cell" buttons), though the prose shows as comments.

### Pip-only alternative (no conda)

All dependencies are on PyPI, so a plain virtual environment works too:

```bash
python -m venv env && source env/bin/activate
pip install "ipi>=3.2.0" numpy matplotlib
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
├── scripts/                    # i-PI launching plumbing (run_ipi/run_parallel) + plotting helpers
└── data/                       # the five i-PI input files (3- and 4-particle cases)
```

The example directory follows the
[atomistic-cookbook](https://github.com/lab-cosmo/atomistic-cookbook)
conventions (a `README.rst`, an `environment.yml`, and a sphinx-gallery `.py`
that renders to a notebook + HTML page).

## The exact benchmarks

Each case is compared against the **exact** energy of non-interacting particles
in a harmonic trap, computed in `analysis.py` from the canonical
partition-function recursion (mHa):

| Case | βℏω₀ | exact energy (mHa) |
|------|------|--------------------|
| 4 distinguishable | 2.00 | 0.869 |
| 4 bosons | 2.00 | 0.749 |
| 3 bosons + 1 dist | 2.00 | 0.798 |
| 3 fermions | 1.16 | **1.0530** |

(The statistics comparison uses four particles; the fermion case uses three. The
boson energy–temperature sweep also uses three bosons, printed as `E / ℏω₀` at
each `βℏω₀`.)

## Fermions need averaging and careful error bars

Fermionic averages come from simulating **bosons** and reweighting each sample by
the fermionic sign `s`: `⟨A⟩_F = ⟨A s⟩ / ⟨s⟩`. When the average sign is small —
the **sign problem** — this ratio has a large variance, so the tutorial runs
several short trajectories and combines them with the **sign-weighted** estimator
from the SI of Hirshberg, Invernizzi & Parrinello (*JCP* **152**, 171102, 2020).
Each trajectory `j` is weighted by its total sign `W_j = Σ s`, with an effective
sample size `n_eff`:

```
Ē_F = Σ Wⱼ Eⱼ / Σ Wⱼ ,   n_eff = (Σ Wⱼ)² / Σ Wⱼ² ,   error = σ_E / √n_eff
```

Trajectories carry unequal information (a small-sign trajectory has a
poorly-determined ratio), so `n_eff < M` and the honest error bar is a little
larger than `std/√M`. `weighted_average()` in `analysis.py` reduces to the plain
mean + `std/√M` when the weights are equal, so the same code handles the
sign-free bosons. At the tutorial's settings (three fermions, 30 K, βℏω₀ = 1.16,
8 trajectories) this gives **1.06 ± 0.05 mHa vs the exact 1.053** (n_eff ≈ 7 of
8) — bracketing the exact value.

The tutorial deliberately runs the fermions warm (30 K) and short. Going colder,
or to more particles, shrinks the sign and blows up the per-trajectory scatter —
the recipe shows this by repeating the run at βℏω₀ = 2, where the eight
trajectories fly apart. For production you would run more (and longer)
trajectories and, if needed, more beads.

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
