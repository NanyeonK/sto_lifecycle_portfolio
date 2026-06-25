#!/usr/bin/env bash
# Run E1_2L baseline under v4 Option 1 state extension.
# Coarse grids: N_W=15, N_Z=5, N_X_PREV=3 (x_prev ∈ {0, 0.5, 1.0}).
# Wall time estimate: ~1.5-2h single thread.
set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

REGIME=E1_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
RHO_AB=0.50 P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e1.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e1_stdout.log"

echo "E1_2L v4 done. Summary at $OUTDIR/p6_option1_e1.json"
