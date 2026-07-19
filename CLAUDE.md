# CLAUDE.md — project state & context

Working notes for continuing this project (readable by Claude Code on any
machine). This repo is a **bosons/fermions PIMD tutorial for i-PI 3.x**,
structured as an [Atomistic Cookbook](https://atomistic-cookbook.org) recipe,
prepared for the **PIQM 2026** meeting in Shanghai. It modernises the 2023
tutorial module [`i-pi/piqm2023-tutorial/06-indistinguishable`](https://github.com/i-pi/piqm2023-tutorial/tree/main/06-indistinguishable).

Owner: Barak Hirshberg. Cookbook conventions modelled on
`lab-cosmo/atomistic-cookbook` (per-example `environment.yml` + sphinx-gallery
`.py` + `data/`).

## Layout

```
examples/bosons-fermions-pimd/
  bosons-fermions-pimd.py   # the tutorial (sphinx-gallery format)
  analysis.py               # i-PI output reader, EXACT analytical energies, fermionic reweighting + SI weighted error + block_average
  environment.yml           # conda: python + pip:[ipi>=3.2, numpy, matplotlib, ase, chemiscope]
  data/*.xml                # 4 i-PI input TEMPLATES: 3dist / 3bosons / 2bosons1dist / 3fermions (recipe overrides T/nbeads/seed/steps per run)
README.md, noxfile.py, LICENSE(BSD-3)
SM.pdf                      # Barak's paper SI (gitignored, NOT redistributed) — source of the fermion error method
```

## Physics / setup

3 particles, mass 1, in a 3D isotropic harmonic trap, k = 1.21647924e-8
Ha/Bohr² (ℏω₀ = 3 meV). Temperatures are quoted as **βℏω₀ = ℏω₀/k_BT**. Four
template inputs in `data/`, each a one-line `<bosons>` change: `[]` dist,
`[0,1,2]` bosons, `[0,1]` mixed, `[0,1,2]` + `fermionic_sign` for fermions
(reweighted). A single `run_ipi()` helper overrides temperature / nbeads / seed /
total_steps per run, so the templates are just starting points.

**Tutorial structure** (redesigned 2026-07-16, Barak's call — three parts):
1. **Bosons + E(T) curve**: sweep βℏω₀ ∈ {1,2,3,5} with *temperature-scaled*
   beads P ∈ {8,16,24,32} (Trotter error ~(βℏω₀/P)², so P grows with βℏω₀ — keeps
   every point ~0.2% converged and the warm point cheap). 4 runs in parallel,
   each with a **block-averaged error bar**; plotted on the exact boson curve +
   dashed ground-state line (bosons/dist → 4.5 ℏω₀, fermions → 6.5 ℏω₀).
2. **Switching statistics** at βℏω₀=2 (17.4 K), P=16: dist/bosons/mixed run in
   parallel, bar chart with error bars vs exact.
3. **Fermions** at βℏω₀=1.16 (30 K), P=12: single reweighted run (noisy) → 8-traj
   sign-weighted ensemble (brackets exact) → SI error-estimation discussion.

Sweep/switch runs use **total_steps=6000, skip=1000** (Barak kept 6000 for a
clean monotonic curve, 2026-07-16); fermion runs use the template 3000 steps.

Driver: pip `i-pi-py_driver -m harmonic -o k` (no compilation). The recipe
**auto-detects** the compiled f90 `i-pi-driver -m harm3d` if on PATH. **NOTE
(measured 2026-07-16): for this tiny 3-atom system the f90 driver is NOT
meaningfully faster** — i-PI's own per-step Python overhead dominates, not the
force eval (~3000 bead-steps/s either way). So "run longer for a cleaner curve"
is expensive regardless of driver; the fix was parallel runs + honest error bars,
NOT more steps. Conda has NO `ipi` package, so f90 is opt-in via compilation
(clone i-pi source, `make -C drivers/f90` → binary at `i-pi-src/bin/i-pi-driver`).

## KEY FINDING (the reason for this whole investigation)

The fermion PIMD "did not agree with the exact result." **Root cause: the
analytical fermion reference was WRONG.** The 2023 code hard-coded a closed form
giving 0.912 mHa at 30 K. The correct exact value is **1.053 mHa**, confirmed by
TWO independent methods (canonical symmetric-polynomial recursion + brute-force
enumeration of 3-fermion states). Fixed in `analysis.py::get_harmonic_energy`
(now the exact recursion for any N, ξ=±1). dist/boson/mixed references were
already correct (0.6514 / 0.5803 / 0.6235 mHa). With the corrected reference,
well-sampled PIMD (~1.05–1.11 mHa) agrees with exact within error — it was a
wrong target, not a simulation bug.

## Fermion error estimation (from SM.pdf, Sec. II)

Implemented in `analysis.py`: per trajectory j, `E_j = <εs>_j/<s>_j` and weight
`W_j = Σ(instantaneous sign)`; weighted mean `Ē_F = ΣW_jE_j/ΣW_j` (= pooled
estimator), effective sample size `n_eff = (ΣW)²/ΣW²`, error `σ/√n_eff`
(Ambegaokar–Troyer). `weighted_average()` reduces to plain mean + std/√M for
equal weights (bosons, average sign). Checked (20 traj × 5000): naive
mean-of-ratios is biased high (1.45) vs weighted 1.11; n_eff-based error is ~12%
larger than assuming n=M. There is a dedicated recipe + README section on this.

## Status (handoff — continue on a faster machine)

Everything below is committed and pushed to
`https://github.com/BarakHirshberg/PIQM2026_BH_Tutorial` (branch `master`).
Working tree is clean; `git pull` to get the latest.

DONE:
- [x] **Root cause found & fixed**: wrong hard-coded 3-fermion benchmark
      (0.912 → exact **1.053 mHa**), verified two independent ways. `analysis.py`
      now uses the exact recursion (ξ=±1) for any N.
- [x] SI weighted fermion error estimation in `analysis.py`
      (`weighted_average`, `fermionic_trajectory_estimate`), unit-checked.
- [x] f90 driver: build + validate + recipe **auto-detection**
      (`driver_command()` uses `i-pi-driver -m harm3d` if on PATH, else
      `i-pi-py_driver -m harmonic`).
- [x] Socket-race fix (poll for `/tmp/ipi_<addr>` before launching the driver,
      which does not retry) in the recipe's `_run_seeded`.
- [x] **Tutorial redesigned** per Barak's decisions (2026-07-16): NOT tightly
      converged; keep it laptop-light and "make sense." Fermion section = one
      short run (flagged unreliable) + a **light parallel multi-trajectory**
      average (`run_fermion_ensemble`, 8 short trajectories, sign-weighted)
      whose large-but-honest error bar brackets the exact value. All 0.912 →
      1.053 fixed in recipe + README.
- [x] **Removed the whole `reference/` directory** (Barak's call 2026-07-16):
      with exact analytical benchmarks and the inline light multi-trajectory
      demo, the heavy well-sampled reference table was redundant. The SI error
      estimator lives in `analysis.py` and is used inline by the recipe.
- [x] **Tutorial RE-STRUCTURED into 3 parts** (2026-07-16, this session): bosons
      + E(T) sweep → switching statistics → fermions (see Physics/setup above).
      Unified `run_ipi()` helper (overrides T/nbeads/seed/steps, isolated tmpdir,
      BLAS-thread cap). Boson sweep + switching run in PARALLEL with block-averaged
      error bars (`analysis.block_average`, `analysis.total_energy_series`). README
      rewritten to match (βℏω₀ units, trap freq, bead-scaling table).
- [x] **VERIFIED end-to-end** (2026-07-16, f90 driver, full run, exit 0): boson
      sweep monotonic onto the 4.5 ℏω₀ ground state; switching bars ordered
      bosons<mixed<dist, ~1σ of exact; **fermion 8-traj weighted 1.064 ± 0.052 mHa
      brackets exact 1.053 within 0.2σ** (n_eff 6.9/8). All 3 figures render.

REMAINING — one thing:
- [ ] **Noob-grad-student test**: spawn an agent role-playing an inexperienced
      i-PI user; have it create a FRESH env (`./env_noobtest`) following ONLY
      README.md, run the whole recipe, and report friction points; then fix them
      and remove the test env. (Fresh env tests the pip-only path honestly.)

Runtime note: sweep/switch at 6000 steps ≈ 5–6 min on a 20-core f90 machine, but
~10–15 min on a 4-core laptop with the pip Python driver (f90 barely helps — see
Physics/setup). Barak accepted this for the clean curve (2026-07-16).

Known open question (minor, not blocking): short runs scatter a few % around
exact; the block-averaged error bars are what make that "make sense." Tightening
would need many more (or longer) steps — expensive, and not the point.

## Session 2026-07-18 — units, propagator, block averaging (all on master, pushed)

All committed & pushed to `master` this session (`git pull` to get them):
- **Energies now in ℏω₀ everywhere** (was mHa in places): all 4 figures
  (`scripts/plots.py`), the switching + fermion **printouts** in the recipe, and
  the README benchmark table (7.88 / 6.79 / 7.23 / **9.55**). Recorded under
  ## Conventions. The one mHa kept parenthetically is the 1.053 mHa fermion
  benchmark = 9.55 ℏω₀. `plots.py` also renamed unused `label`→`_label`.
- **README cleanups**: "(mass 1)" → "$m=1$"; corrected the f90 "a few times
  faster" claim (it is NOT faster for this tiny system); pruned 2 dead inputs
  (`input_3dist`, `input_2bosons1dist`) → `data/` now has 5 files ("seven"→"five").
- **Propagator note added** (tutorial intro): bosons REQUIRE `propagator='bab'`.
  Verified in i-PI **3.3.0**: `engine/normalmodes.py::free_qstep` raises
  `NotImplementedError` for `exact`/`cayley` when `<bosons>` is set, and the input
  default propagator is `exact`, so bab must be set explicitly (doc:
  `inputs/normalmodes.py` propagator help).
- **Atom ordering is NOT constrained** in i-PI 3.3.0 (checked; no change needed):
  `utils/exchange.py` reads boson coords via fancy indexing
  `beads.q.reshape(...)[:, self.bosons]`, masses `beads.m[self.bosons]`, forces
  written back at `fspring[:, self.bosons, :]`. Bosons may be any indices,
  contiguous or not; distinguishables may precede them.
- **block_average n_blocks=10 justified** (was arbitrary): Flyvbjerg–Petersen
  blocking on fresh trajectories → total-energy autocorrelation time ~25–55 steps,
  so ~500-step blocks (10 blocks of ~5000 samples) ≈ 10–20 τ. Documented in the
  docstring; kept FIXED (not recomputed per call, per Barak). Warmest points
  (βℏω₀=1,2, incl. the switching condition) are near the resolution limit of a
  6000-step run — their bars are mild lower bounds.
- **Content review done** (two agents): the physics reviewer VERIFIED the tutorial
  is correct (4.5/6.5 ℏω₀ ground states, 1.053 mHa, conversions, ordering,
  estimator math). Pedagogy findings were deferred (Barak: "not serious"). The
  full list is parked in the plan file `humble-napping-wall.md`. Don't redo.

**Notebook viewer (was running this session in /tmp — GONE after a restart).**
A fully-executed `bosons-fermions-pimd.ipynb` was served via `jupyter lab` on
`127.0.0.1:8888` (token `piqm2026`) and viewed over SSH/Termius port-forwarding.
The server and its conda env live under the session scratchpad (`.nox/…`), which
`/tmp` wipes on restart. To view again in a fresh session: get an env with
`i-pi + numpy + matplotlib + jupyterlab + sphinx-gallery`, then
`sphinx_gallery_py2jupyter bosons-fermions-pimd.py`,
`jupyter nbconvert --to notebook --execute --inplace bosons-fermions-pimd.ipynb`
(~15 min), then `jupyter lab --no-browser --port=8888 --ip=127.0.0.1
--ServerApp.token=piqm2026 --ServerApp.root_dir=<example dir>` and port-forward
8888 in Termius. (The `.ipynb` is gitignored; it is a build artifact.)

## Cookbook contribution (Phase 2) — #292 MERGED, follow-up #294 open

Submitted upstream as **[lab-cosmo/atomistic-cookbook#292](https://github.com/lab-cosmo/atomistic-cookbook/pull/292)**
(base `main`, head `BarakHirshberg:add-bosons-fermions-pimd`, 14 files, +1391/-0).
**Validated end-to-end** with the real harness on a fast machine
(`nox -e bosons-fermions-pimd` via miniconda's conda 25.x): the example session +
`build_website` both succeed, recipe runs in ~2–3 min, no example-specific sphinx
warnings. CI on the PR waits for a maintainer "approve and run" (first-time fork
contributor).

**#292 MERGED (2026-07-18, merge commit `6e99c2f`).** Its lint job failed first
(the cookbook runs flake8 `--max-line-length=88` + flake8-bugbear + blackdoc +
`ruff format --check` — stricter than `ruff check`: 3× E501, 1× B007); Michele
Ceriotti (`ceriottm`) pushed a fix ("Please the lint gods", `204ce238`) that
introduced a `centroid-virial` typo, Barak pushed the one-word fix, and it merged.
**During review the helper modules were moved into `data/`** (not `scripts/`):
the merged example has `data/analysis.py`, `data/plots.py`, `data/ipi_runs.py`,
`data/__init__.py`, imported as `from data import analysis, plots` /
`from data.ipi_runs import run_parallel` (`plots.py` uses `from . import analysis`).
Lint like CI locally with `flake8 --max-line-length=88 --extend-ignore=E203` +
`ruff format --check` + `blackdoc --check` (NOT just `ruff check`).

**Master's 2026-07-18 changes were synced into the merged cookbook via follow-up
[PR #294](https://github.com/lab-cosmo/atomistic-cookbook/pull/294)** (branch
`sync-hw0-units-propagator` off upstream/main, 3 files, +42/-20, lint-clean):
ℏω₀ units in the remaining figures + printouts, the `propagator='bab'` note, and
the block_average one-line note. Cookbook-only deltas (1500/2000 steps, the
step-count note, `data/` layout, trimmed inputs) preserved. OPEN, awaiting CI/merge.
Fork clone with both branches lives at `~/atomistic-cookbook` (origin=fork,
upstream=lab-cosmo).

Cookbook version diverges from master on purpose (master kept untouched):
- shorter steps for the 12-min CI budget (`SWEEP_STEPS=1500`, `FERMION_STEPS=2000`,
  `SKIP=300`) + a `.. note::` telling readers to increase them;
- dropped `ase`/`chemiscope` from `environment.yml` (unused);
- pruned 2 dead 3-particle inputs (`input_3dist`, `input_2bosons1dist`);
- **helper modules live in `data/`** (`data/analysis.py`, `data/plots.py`,
  `data/ipi_runs.py`, `data/__init__.py`) — a maintainer moved them there during
  the #292 review so sphinx-gallery (pattern `.*`, top-level `.py` only) renders
  just the one recipe `.py`. Imports: `from data import analysis, plots` /
  `from data.ipi_runs import run_parallel` (`plots.py`: `from . import analysis`).
  Master keeps `analysis.py` top-level + `scripts/` for the rest; layouts differ.

Fork/PR mechanics on THIS machine: nox needs Python ≥3.9 (`list[Path]` in
`src/get_examples.py`) — build the nox venv with `/home/hirshb/miniconda/bin/python`
(3.13), NOT the 3.8 conda envs. `venv_backend="conda"` needs a conda supporting
`--solver=libmamba` (miniconda's 25.x, not anaconda3's old conda) with accepted
channel ToS (`conda tos accept` for pkgs/main + pkgs/r, already done). Persistent
fork clone at `~/atomistic-cookbook`; the noxvenv for local lint is under the
session scratchpad (wiped on restart).

**Resuming cookbook work later** (#292 is merged; the example lives on upstream
`main`). For further changes, branch off **upstream/main** and open a **new** PR
(as #294 did) — do NOT push to the merged `add-bosons-fermions-pimd`:
```bash
cd ~/atomistic-cookbook          # origin=fork, upstream=lab-cosmo
git fetch upstream main && git checkout -b <topic> upstream/main
# edit under examples/bosons-fermions-pimd/ (helpers in data/), commit, then:
git push origin <topic>
gh pr create --repo lab-cosmo/atomistic-cookbook --base main \
  --head BarakHirshberg:<topic>
```
Lint like CI before pushing (`flake8 --max-line-length=88 --extend-ignore=E203`
+ `ruff format --check` + `blackdoc --check`). To re-validate the build: nox venv
with `/home/hirshb/miniconda/bin/python`, miniconda's conda on PATH,
`nox -e bosons-fermions-pimd`.

## Gotchas / lessons
- Soft trap: oscillation period ~1370 fs, so dt=1 fs needs millions of steps.
  Internal ring-polymer modes are propagated exactly (normal modes + nmts=10),
  so **dt=10 fs** is safe and ~10× cheaper. All inputs use dt=10.
- The SI weighted mean with `W_j = Σs` equals the pooled `Σεs/Σs` estimator.
- f90 driver does NOT retry the socket connection; under concurrency you must
  wait for `/tmp/ipi_<address>` to exist before launching it (done in `run_ipi`).
- For this tiny 3-atom system the f90 driver is ~as slow as the Python one —
  i-PI per-step overhead dominates (~3000 bead-steps/s). Don't reach for longer
  runs to clean up noise; parallelise + use block-averaged error bars instead.
- /tmp (incl. the scratchpad f90 build) is cleared on session restart; rebuild
  with: `git clone --depth 1 https://github.com/i-pi/i-pi && make -C i-pi/drivers/f90`
  (binary lands at `i-pi/bin/i-pi-driver`; put its dir on PATH).

## Conventions
- **Energy units: always ℏω₀** in every figure AND every printout (divide the Ha
  value by `analysis.omega0()`), not mHa — so all panels/outputs are directly
  comparable. The plot helpers in `scripts/plots.py` convert internally; the
  recipe's `print()`s do it inline. The README benchmark table is in ℏω₀ too. The
  one mHa value kept (parenthetically) is the 1.053 mHa fermion benchmark
  (= 9.55 ℏω₀), since it's the recognizable KEY-FINDING number.
- Verify with the conda env at `./env` (gitignored). Activate:
  `source ~/miniconda3/etc/profile.d/conda.sh && conda activate ./env`.
  (On a machine without that conda, make a scratch venv `./env_verify` —
  gitignored — `python -m venv env_verify && env_verify/bin/pip install
  "ipi>=3.2" numpy matplotlib ruff`; that's what this session used.)
- Lint: `ruff check` + `ruff format --check`. All committed code passes both.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
- Do NOT commit SM.pdf (gitignored) or the ./env folder.
