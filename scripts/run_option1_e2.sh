#!/usr/bin/env bash
# Run E2_2L baseline with v4 (Option 1 full state extension).
# Default grids: N_W=15, N_Z=5, N_X_PREV=3.
# Approx wall: ~70-120 min single thread.
# Run from repo root: bash scripts/run_option1_e2.sh

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

REGIME=E2_2L \
N_W=${N_W:-15} N_Z=${N_Z:-5} \
N_X_PREV=${N_X_PREV:-3} X_PREV_MAX=${X_PREV_MAX:-1.5} \
ASSET_GRID_SIZE=${ASSET_GRID_SIZE:-9} X_GRID_SIZE=${X_GRID_SIZE:-5} \
GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo "E2_2L v4 done. Summary: $OUTDIR/p6_option1_e2.json"
