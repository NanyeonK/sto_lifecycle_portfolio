#!/usr/bin/env bash
# Run E2_2L (Option 1 full state extension) baseline at v4 grid settings.
# Execute on server1 in the sto_lifecycle_portfolio working directory.
# Usage:  bash scripts/run_option1_e2.sh
#
# Produces: output/diagnostics/p6_option1_e2.json + stdout log
# Expected wall time: ~2-2.5 hours (single thread, N_W=15 N_Z=5 N_X_PREV=3)

set -e
cd "$(dirname "$0")/.."

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "[run_option1_e2] Starting E2_2L v4 at $(date)"

REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo "[run_option1_e2] Done at $(date)"
