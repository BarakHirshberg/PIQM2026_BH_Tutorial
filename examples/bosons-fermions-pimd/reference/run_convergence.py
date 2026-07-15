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
CASE_FILTER = os.environ.get("CASE_FILTER")  # substring to select a subset of cases
SPRING_CONSTANT = "1.21647924e-8"
SKIP = int(os.environ.get("SKIP", "200"))        # thermalisation rows to discard
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
        drv = subprocess.Popen(
            ["i-pi-py_driver", "-m", "harmonic", "-o", SPRING_CONSTANT, "-u", "-a", address],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        ipi.wait()
        drv.wait()

        out = os.path.join(tmp, "data.out")
        if is_fermionic:
            sign, energy = analysis.reweighted_fermionic_energy(out, SKIP)
            return energy, sign
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
            energies = np.array([r[0] for r in results])
            mean = energies.mean()
            sem = energies.std(ddof=1) / np.sqrt(len(energies))
            ref = analysis.analytical_energy(temp, sys_type)
            extra = ""
            if is_fer:
                signs = np.array([r[1] for r in results])
                extra = f"  <sign> = {signs.mean():.3f} +/- {signs.std(ddof=1)/np.sqrt(len(signs)):.3f}"
            rows.append((label, mean, sem, ref))
            print(
                f"{label:20s} {mean*1e3:7.4f} +/- {sem*1e3:6.4f} mHa   "
                f"(analytical {ref*1e3:7.4f} mHa){extra}   "
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
