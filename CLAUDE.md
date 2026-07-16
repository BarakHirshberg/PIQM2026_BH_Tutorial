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
  analysis.py               # i-PI output reader, EXACT analytical energies, fermionic reweighting + SI weighted error
  environment.yml           # conda: python + pip:[ipi>=3.2, numpy, matplotlib, ase, chemiscope]
  data/*.xml                # 4 i-PI inputs: 3dist / 3bosons / 2bosons1dist / 3fermions
README.md, noxfile.py, LICENSE(BSD-3)
SM.pdf                      # Barak's paper SI (gitignored, NOT redistributed) — source of the fermion error method
```

## Physics / setup

3 particles, mass 1, in a 3D isotropic harmonic trap, k = 1.21647924e-8
Ha/Bohr² (ℏω₀ = 3 meV). Four cases, each a one-line `<bosons>` change:
`[]` dist, `[0,1,2]` bosons, `[0,1]` mixed (2 bosons + 1 dist), `[0,1,2]` +
`fermionic_sign` for fermions (reweighted). dist/bosons/mixed at 17.4 K, P=32;
fermions at 30 K, P=12 (higher T / fewer beads because of the sign problem).

Driver: pip `i-pi-py_driver -m harmonic -o k` (no compilation). The recipe
**auto-detects** the compiled f90 `i-pi-driver -m harm3d` if on PATH (~3×
faster; identical results). Conda has NO `ipi` package, so the f90 driver is
opt-in via compilation (clone i-pi source, `make -C drivers/f90`).

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

REMAINING — exactly two things, do on the faster machine:
- [ ] **1. VERIFY the recipe end-to-end** — NOT yet confirmed (the run was cut
      off here). Run it and check the new multi-trajectory cell works and the
      numbers make sense:
      `cd examples/bosons-fermions-pimd && MPLBACKEND=Agg python bosons-fermions-pimd.py`
      (put the f90 `i-pi-driver` on PATH to make it fast). Expected: dist/boson/
      mixed single runs near the exact values; fermion 8-traj weighted energy
      ~1.0–1.1 mHa with a large error bar that brackets 1.053; two figures render.
- [ ] **2. Noob-grad-student test**: spawn an agent role-playing an
      inexperienced i-PI user; have it create a FRESH env (`./env_noobtest`)
      following ONLY README.md, run the whole recipe, and report friction points;
      then fix them and remove the test env. (Fresh env tests the pip-only path
      honestly.)

Known open question (minor, not blocking): with tight statistics all cases sit
~1–2% (≈2σ) above exact — a small finite-P(32)/dt(10 fs) systematic. Barak is
fine with this; the tutorial is meant to "make sense," not be converged. If you
ever want <1σ agreement, use more beads (P≥64) and/or smaller dt in the
reference runs only.

## Gotchas / lessons
- Soft trap: oscillation period ~1370 fs, so dt=1 fs needs millions of steps.
  Internal ring-polymer modes are propagated exactly (normal modes + nmts=10),
  so **dt=10 fs** is safe and ~10× cheaper. All inputs use dt=10.
- The SI weighted mean with `W_j = Σs` equals the pooled `Σεs/Σs` estimator.
- f90 driver does NOT retry the socket connection; under concurrency you must
  wait for `/tmp/ipi_<address>` to exist before launching it.
- /tmp (incl. the scratchpad f90 build) is cleared on session restart; rebuild
  with: `git clone --depth 1 https://github.com/i-pi/i-pi && make -C i-pi/drivers/f90`.

## Conventions
- Verify with the conda env at `./env` (gitignored). Activate:
  `source ~/miniconda3/etc/profile.d/conda.sh && conda activate ./env`.
- Lint: `ruff check`. All committed code passes ruff.
- Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
- Do NOT commit SM.pdf (gitignored) or the ./env folder.
