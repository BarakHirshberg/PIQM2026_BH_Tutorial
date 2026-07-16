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
  reference/
    run_convergence.py      # multi-trajectory study: mean +/- error, SI weighting for fermions
    run_final_table.sh      # generates the well-sampled reference table (fermion P-scan)
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

## Status of the "5.5-hour autonomous" task

DONE:
- [x] Found & fixed the wrong analytical reference (committed).
- [x] Built + validated the f90 driver; wired auto-detection into the recipe;
      `IPI_DRIVER` env in run_convergence.py (committed).
- [x] Fixed a hang: run_convergence now polls for the i-PI socket before
      launching the (non-retrying) f90 driver (committed).
- [x] Bead sweep P=12/24/36/48 (light, 10 traj × 8000): values scatter
      1.00–1.09 around exact 1.053, no clean Trotter trend — likely just
      statistical; needs heavier sampling to be clean.

IN PROGRESS / REMAINING (a background job may still be running here):
- [ ] **Final well-sampled table** — `reference/run_final_table.sh` runs the 3
      easy cases (P=32, 24 traj × 8000) + a fermion P-scan (12/24/48, 24 traj ×
      16000) with the f90 driver. Output goes to `/tmp/final_table.log`. When it
      finishes, put the converged numbers into the recipe's "Convergence"
      table (currently a placeholder) and the "Error estimation" section
      (currently still cites the OLD 0.912 exact — MUST update to 1.053), and
      into the README table. Confirm all four cases agree with exact within error.
- [ ] **Tutorial narrative pass**: update every remaining mention of 0.912 →
      1.053; the fermion conclusion changes from "does not overlap / gap is
      physical" to "agrees once the reference is fixed and errors are weighted".
      Lines to fix in bosons-fermions-pimd.py: ~262 (table row) and ~305–322
      (error-estimation numbers/conclusion). Same in README.md fermion section.
- [ ] **End-to-end run** of the recipe with corrected refs + driver auto-detect
      (~1.5 min with f90, ~4.5 min pip driver) to confirm it still works.
- [ ] **Noob-grad-student agent test**: spawn an agent role-playing an
      inexperienced i-PI user; have it create a FRESH conda env from the README
      and run the whole tutorial, reporting friction points; then fix them.
      Draft brief saved at (scratch, may be gone after restart) — re-derive from
      this file: fresh `./env_noobtest`, follow only README, report + clean up.

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
