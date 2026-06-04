#!/usr/bin/env bash
# Run E1_2L baseline with v4 (6D state, Option 1 tx_cost on deltas)
# Default coarse grids: N_W=15, N_Z=5, N_X_PREV=3
# Expected wall time: ~2-4 hours per regime on server1 (single thread)

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "=== v4 E1_2L baseline run ==="
date

REGIME=E1_2L \
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
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e1.json" \
JULIA_NUM_THREADS=1 \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e1_stdout.log"

echo "=== E1_2L done ==="
date
