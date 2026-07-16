#!/usr/bin/env bash
# Generate the final, well-sampled reference table for all four cases.
#
# The distinguishable/bosonic cases converge quickly; the fermionic case needs
# far more sampling (the sign problem) and is run as a P-scan (12/24/48 beads)
# to show it converges to the exact value. Uses the compiled f90 driver if
# IPI_DRIVER points to it (a few times faster), otherwise the pip driver.
#
#   IPI_DRIVER=/path/to/i-pi-driver bash reference/run_final_table.sh
set -u
cd "$(dirname "$0")/.."
LOG=${LOG:-/tmp/final_table.log}
: > "$LOG"

run() {  # label_filter  NBEADS  STEPS  NTRAJ  SKIP
  echo "===== $1 (P=$2, $4 traj x $3 steps) =====" | tee -a "$LOG"
  N_TRAJ="$4" STEPS="$3" NBEADS="$2" SKIP="$5" MAX_CONCURRENT="${MAX_CONCURRENT:-12}" \
    CASE_FILTER="$1" python reference/run_convergence.py 2>&1 | tee -a "$LOG"
}

# distinguishable / bosonic cases: converge fast, 24 traj x 8000 steps @ P=32
run "distinguishable" 32 8000 24 400
run "3 bosons"        32 8000 24 400
run "1 dist"          32 8000 24 400
# fermions: sign problem -> heavy sampling; P-scan to show convergence to exact
run "fermion"         12 16000 24 600
run "fermion"         24 16000 24 600
run "fermion"         48 16000 24 600

echo "FINAL TABLE DONE" | tee -a "$LOG"
