#!/bin/bash
# run_option1_e2.sh — E2_2L baseline for v4 Option 1 state extension
# Proper tau_buy via x_prev state: pre-holding x_B at ell=A saves tau_buy at relocation.
# Run on server1 from repo root:
#   bash scripts/run_option1_e2.sh
#
# Output: output/diagnostics/p6_option1_e2.json

set -euo pipefail

OUTDIR=output/diagnostics
mkdir -p "$OUTDIR"

echo "=== v4 E2_2L Option 1 baseline ==="
echo "Start: $(date)"

REGIME=E2_2L \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
JULIA_NUM_THREADS=1 \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo "End: $(date)"
echo "Summary: $OUTDIR/p6_option1_e2.json"
