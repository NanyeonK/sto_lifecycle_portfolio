#!/usr/bin/env bash
# Run E1_2L baseline with v4 6D-state solver (Option 1).
# Grids: N_W=15, N_Z=5, N_X_PREV=3 (coarse; ~2-4h wall on server1).
# Compare V_t1_midpoint_ellA_xprev0 with v3 E1_2L result to confirm consistency.
# Usage: bash scripts/run_option1_e1.sh

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

REGIME=E1_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=5 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e1.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e1_stdout.log"

echo "E1_2L v4 done. Summary: $OUTDIR/p6_option1_e1.json"
