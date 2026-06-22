#!/usr/bin/env bash
# Run E2_2L baseline with v4 (6D state, Option 1 tx_cost on deltas)

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "=== v4 E2_2L baseline run ==="
date

REGIME=E2_2L \
N_W="${N_W:-15}" \
N_Z="${N_Z:-5}" \
N_X_PREV="${N_X_PREV:-3}" \
X_PREV_MAX="${X_PREV_MAX:-1.5}" \
ASSET_GRID_SIZE="${ASSET_GRID_SIZE:-7}" \
X_NEW_GRID_SIZE="${X_NEW_GRID_SIZE:-5}" \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
SAVE_PATH="$OUTDIR/p6_option1_e2_result.jls" \
JULIA_NUM_THREADS=1 \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo "=== E2_2L done ==="
date
