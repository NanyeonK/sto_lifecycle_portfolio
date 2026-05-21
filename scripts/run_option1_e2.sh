#!/usr/bin/env bash
# run_option1_e2.sh — E2_2L baseline run with v4 Option 1 state extension
# Usage: bash scripts/run_option1_e2.sh
# Prereq: Julia installed, JSON3 available in global Julia env.
#
# Key output: mean_xB_t1_init_ellA — if > 0, hedge channel activated.
# Compare with v3 E2_2L result (mean_xB was 0 there) to verify the
# delta tau_buy mechanism resurrects cross-location holding.

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "[$(date)] Starting E2_2L v4 Option 1 run"

REGIME=E2_2L \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
N_W=15 \
N_Z=5 \
ASSET_GRID_SIZE=9 \
GH_NODES=3 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.005 \
APPLY_TAU_BUY=1 \
RHO_AB=0.50 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
  julia src/vfi_solver_v4.jl \
  2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo "[$(date)] E2_2L v4 done. Output: $OUTDIR/p6_option1_e2.json"
echo "Check: mean_xB_t1_init_ellA in p6_option1_e2.json should be > 0 if hedge channel active"
