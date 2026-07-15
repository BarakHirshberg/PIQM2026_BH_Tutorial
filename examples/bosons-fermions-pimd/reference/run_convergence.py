"""Offline convergence study: many independent trajectories per case.

A single short i-PI run (as in the tutorial) is noisy. Because the dynamics are
NVT with a proper thermostat there is *no systematic bias* -- averaging several
independent trajectories (different random seeds) converges towards the exact
value, and the spread across trajectories gives a proper statistical error.

This script runs ``N_TRAJ`` trajectories for each of the four cases, using the
*same* short settings as the tutorial (so the point is purely "average over
trajectories"), parallelised across CPU cores, and prints a Markdown table of

    mean total energy  +/-  standard error of the mean   vs   analytical.

It is meant to be run offline (it takes a few minutes); the resulting table is
copied into the recipe. Run from the example directory::

    python reference/run_convergence.py

Nothing here changes the pip-only tutorial itself -- it only uses the same
``i-pi`` + ``i-pi-py_driver`` already installed.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

# make the sibling analysis.py importable when run from the example directory
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import analysis  # noqa: E402


# --------------------------------------------------------------------------- #
#  Configuration                                                              #
# --------------------------------------------------------------------------- #
# All of these can be overridden from the environment, e.g.
#   N_TRAJ=20 STEPS=10000 CASE_FILTER=fermion python reference/run_convergence.py
N_TRAJ = int(os.environ.get("N_TRAJ", "10"))     # independent trajectories per case
MAX_CONCURRENT = int(os.environ.get("MAX_CONCURRENT", "8"))  # simultaneous i-PI+driver pairs (<= cores)
STEPS = os.environ.get("STEPS")   # override total_steps for every case (None = use the XML value)
NBEADS = os.environ.get("NBEADS")  # override nbeads for every case (None = use the XML value)
CASE_FILTER = os.environ.get("CASE_FILTER")  # substring to select a subset of cases
SPRING_CONSTANT = "1.21647924e-8"
SKIP = int(os.environ.get("SKIP", "200"))        # thermalisation rows to discard
# Force driver: by default the pip-installed pure-Python driver. Set
# IPI_DRIVER to the path of the compiled f90 `i-pi-driver` (mode harm3d) for a
# ~3x speed-up when generating reference data offline.
F90_DRIVER = os.environ.get("IPI_DRIVER")
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# (tag, input file, analytical type, temperature, is_fermionic)
CASES = [
    ("3 distinguishable", "input_3dist.xml", "dist", 17.4, False),
    ("3 bosons", "input_3bosons.xml", "bosonic", 17.4, False),
    ("2 bosons + 1 dist", "input_2bosons1dist.xml", "mixed", 17.4, False),
    ("3 fermions", "input_3fermions.xml", "fermionic", 30.0, True),
]


def run_one(xml_name, seed, is_fermionic):
    """Run a single trajectory with a given seed in an isolated temp dir.

    Returns the total energy (Ha); for the fermionic case this is already the
    sign-reweighted energy, and the trajectory's average sign is returned too.
    """
    tmp = tempfile.mkdtemp(prefix="ipi_conv_")
    address = f"conv-{seed}"
    try:
        # template the input: unique seed + unique socket address
        with open(os.path.join(DATA_DIR, xml_name)) as f:
            xml = f.read()
        import re

        xml = re.sub(r"<seed>.*?</seed>", f"<seed> {seed} </seed>", xml, flags=re.S)
        xml = re.sub(r"<address>.*?</address>", f"<address>{address}</address>", xml, flags=re.S)
        if STEPS:
            xml = re.sub(r"<total_steps>.*?</total_steps>", f"<total_steps> {STEPS} </total_steps>", xml, flags=re.S)
        if NBEADS:
            xml = re.sub(r'nbeads="\d+"', f'nbeads="{NBEADS}"', xml)
        xml_path = os.path.join(tmp, "input.xml")
        with open(xml_path, "w") as f:
            f.write(xml)

        sock = f"/tmp/ipi_{address}"
        if os.path.exists(sock):
            os.remove(sock)

        ipi = subprocess.Popen(
            ["i-pi", "input.xml"], cwd=tmp,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        time.sleep(2)
        if F90_DRIVER:
            drv_cmd = [F90_DRIVER, "-m", "harm3d", "-o", SPRING_CONSTANT, "-u", "-a", address]
        else:
            drv_cmd = ["i-pi-py_driver", "-m", "harmonic", "-o", SPRING_CONSTANT, "-u", "-a", address]
        drv = subprocess.Popen(
            drv_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        ipi.wait()
        drv.wait()

        out = os.path.join(tmp, "data.out")
        if is_fermionic:
            # per-trajectory E_j, W_j (sum of signs), mean sign
            return analysis.fermionic_trajectory_estimate(out, SKIP)
        else:
            return analysis.mean_energies(out, SKIP)["total"], None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    seeds = [10007 + 101 * i for i in range(N_TRAJ)]  # distinct, fixed for reproducibility
    rows = []

    cases = [c for c in CASES if not CASE_FILTER or CASE_FILTER.lower() in c[0].lower()]

    with ThreadPoolExecutor(max_workers=MAX_CONCURRENT) as pool:
        for label, xml, sys_type, temp, is_fer in cases:
            t0 = time.time()
            results = list(pool.map(lambda s: run_one(xml, s, is_fer), seeds))
            ref = analysis.analytical_energy(temp, sys_type)

            if is_fer:
                # results = list of (E_j, W_j, mean_sign_j)
                E = np.array([r[0] for r in results])
                W = np.array([r[1] for r in results])
                signs = np.array([r[2] for r in results])

                M = len(E)
                # naive: unweighted mean of per-trajectory ratios, std/sqrt(M)
                naive_mean, naive_err, _ = analysis.weighted_average(E, None)
                # SI: weight each trajectory by its total sign W_j (Eqs. 4-6)
                w_mean, w_err, n_eff = analysis.weighted_average(E, W)
                # apples-to-apples: same weighted mean, but pretend all M count
                # equally (error = sigma/sqrt(M)) vs the proper sigma/sqrt(n_eff)
                sigma = w_err * np.sqrt(n_eff)
                err_assume_M = sigma / np.sqrt(M)
                # average sign is itself an unweighted average over trajectories
                s_mean, s_err, _ = analysis.weighted_average(signs, None)

                rows.append((label, w_mean, w_err, ref))
                print(
                    f"{label:20s}   [M={M} traj, {time.time()-t0:.0f}s]\n"
                    f"    naive mean-of-ratios : {naive_mean*1e3:7.4f} +/- {naive_err*1e3:6.4f} mHa "
                    f"(biased by low-weight trajectories)\n"
                    f"    SI weighted mean     : {w_mean*1e3:7.4f} mHa   (Eq. 4)\n"
                    f"      error assuming n=M : +/- {err_assume_M*1e3:6.4f} mHa\n"
                    f"      error using n_eff  : +/- {w_err*1e3:6.4f} mHa   "
                    f"(n_eff = {n_eff:.1f} of {M}; larger by x{np.sqrt(M/n_eff):.2f})\n"
                    f"    <sign> = {s_mean:.3f} +/- {s_err:.3f}    analytical {ref*1e3:.4f} mHa"
                )
            else:
                energies = np.array([r[0] for r in results])
                mean, err, _ = analysis.weighted_average(energies, None)
                rows.append((label, mean, err, ref))
                print(
                    f"{label:20s} {mean*1e3:7.4f} +/- {err*1e3:6.4f} mHa   "
                    f"(analytical {ref*1e3:7.4f} mHa)   "
                    f"[{len(energies)} traj, {time.time()-t0:.0f}s]"
                )

    # Markdown table for pasting into the recipe / README
    print("\n### Markdown table\n")
    print(f"| Case | mean of {N_TRAJ} trajectories (mHa) | analytical (mHa) |")
    print("|------|------------------------------|------------------|")
    for label, mean, sem, ref in rows:
        print(f"| {label} | {mean*1e3:.4f} ± {sem*1e3:.4f} | {ref*1e3:.4f} |")


if __name__ == "__main__":
    main()
