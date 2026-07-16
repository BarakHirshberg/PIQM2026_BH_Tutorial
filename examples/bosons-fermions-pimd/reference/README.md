# Reference / production-accuracy studies

These scripts are **not needed to run the tutorial**. They reproduce the
well-sampled numbers behind the tutorial's claims: that once you use the correct
benchmark and average enough trajectories, all four cases agree with the exact
energies, and that the fermionic result converges with the number of beads.

For these heavy runs the compiled Fortran driver is worth building (a few times
faster than the pip Python driver):

```bash
# build once, from a clone of the i-PI source
git clone --depth 1 https://github.com/i-pi/i-pi
make -C i-pi/drivers/f90            # produces i-pi/drivers/f90/driver.x
export IPI_DRIVER=$PWD/i-pi/drivers/f90/driver.x
```

## Scripts

* **`run_convergence.py`** — runs `N_TRAJ` independent trajectories per case in
  parallel (isolated temp dirs, distinct seeds), and reports the mean and a
  proper statistical error. For fermions it uses the sign-weighted estimator
  (Eqs. 4–6 of the SI of Hirshberg, Invernizzi & Parrinello, *JCP* 152, 171102
  (2020)) with an effective sample size. Everything is env-configurable:

  ```bash
  N_TRAJ=20 STEPS=8000 NBEADS=32 MAX_CONCURRENT=6 \
    CASE_FILTER="fermion" IPI_DRIVER=$IPI_DRIVER python run_convergence.py
  ```

* **`run_final_table.sh`** — drives `run_convergence.py` over all four cases with
  well-sampled settings, and does a fermion **bead-convergence scan** (P = 12,
  24, 48).

  ```bash
  IPI_DRIVER=$IPI_DRIVER bash run_final_table.sh   # writes /tmp/final_table.log
  ```

## The exact benchmarks

Computed in `analysis.py` from the canonical partition-function recursion, and
cross-checked by brute-force enumeration of the many-particle states (mHa):

| Case | exact |
|------|-------|
| 3 distinguishable (17.4 K) | 0.6514 |
| 3 bosons (17.4 K) | 0.5803 |
| 2 bosons + 1 dist (17.4 K) | 0.6235 |
| 3 fermions (30 K) | 1.0530 |

The three-fermion value **1.0530** replaces the incorrect 0.912 mHa hard-coded
in the 2023 tutorial — the single biggest reason fermions previously looked like
they disagreed with theory.

## Fermion bead-convergence scan

Sign-weighted energy vs number of beads P (fermions, 30 K; 10 trajectories ×
8000 steps each; exact = 1.053 mHa). The value scatters around the exact result
with no strong bead dependence at this temperature — 12 beads already "make
sense," and the residual scatter is statistical (the sign problem), not a
Trotter trend:

| P (beads) | sign-weighted energy (mHa) | ⟨sign⟩ |
|-----------|----------------------------|--------|
| 12 | 1.09 ± 0.04 | 0.36 |
| 24 | 1.00 ± 0.01 | 0.44 |
| 36 | 1.08 ± 0.05 | 0.34 |

(Re-run with `run_final_table.sh` for tighter, higher-statistics numbers; the
values there use 24 trajectories × 16 000 steps.)

## Error-estimation check

Post-processing 20 trajectories × 5000 steps of the fermion case (exact 1.053
mHa) shows why the weighted estimator matters:

| estimator | mean (mHa) | error (mHa) |
|-----------|-----------|-------------|
| naive mean-of-ratios | 1.45 (biased high) | ± 0.32 |
| SI weighted, assuming n = M | 1.11 | ± 0.078 |
| SI weighted, using n_eff | 1.11 | ± 0.087 (n_eff = 16 of 20) |

The weighting removes the bias in the mean; using the effective sample size
(n_eff = 16 of 20) gives an honest error bar ~12% larger than pretending all 20
trajectories count fully. The weighted result agrees with the exact value.
