#!/usr/bin/env bash
# Run E2_2L baseline under v4 (Option 1: full state extension)
# Coarse grids: N_W=15, N_Z=5, N_X_PREV=3 — intended for first-pass CEV check.
# Expected wall time: ~90-240 min single thread (9x x_prev state factor).
# Run on server1: bash scripts/run_option1_e2.sh

set -euo pipefail
cd "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"

OUTDIR=output/diagnostics
mkdir -p "$OUTDIR"

JULIA_NUM_THREADS=1 \
REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
RHO_AB=0.50 P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
GAMMA=5.0 BETA=0.96 RF=1.02 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo "E2_2L v4 baseline complete. Output: $OUTDIR/p6_option1_e2.json"
