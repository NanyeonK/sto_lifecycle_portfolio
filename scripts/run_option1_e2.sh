#!/usr/bin/env bash
# Run E2_2L baseline under v4 (Option 1 full state extension)
# Coarse grids as per spec: N_W=15, N_Z=5, N_X_PREV=3
set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo "E2_2L v4 done. JSON at $OUTDIR/p6_option1_e2.json"
