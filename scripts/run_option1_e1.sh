#!/usr/bin/env bash
# Run E1_2L baseline at v4 (Option 1) grid settings.
# Execute on server1 in the sto_lifecycle_portfolio working directory.
# Usage:  bash scripts/run_option1_e1.sh
#
# Produces: output/diagnostics/p6_option1_e1.json + stdout log

set -e
cd "$(dirname "$0")/.."

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "[run_option1_e1] Starting E1_2L v4 at $(date)"

REGIME=E1_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e1.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e1_stdout.log"

echo "[run_option1_e1] Done at $(date)"
