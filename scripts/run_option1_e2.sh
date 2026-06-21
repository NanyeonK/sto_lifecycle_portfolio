#!/usr/bin/env bash
# Run E2_2L baseline under v4 (6D state, tau_buy on position deltas)
# Usage: bash scripts/run_option1_e2.sh
# Server1 expected wall time: ~2-3 hours (single thread, N_X_PREV=3)

set -euo pipefail

OUT_DIR="output/diagnostics"
mkdir -p "$OUT_DIR"

REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=9 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 GAMMA=5.0 BETA=0.96 \
SUMMARY_JSON_PATH="$OUT_DIR/p6_option1_e2.json" \
JULIA_NUM_THREADS=1 \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUT_DIR/p6_option1_e2_stdout.log"

echo "E2_2L v4 done. Summary at $OUT_DIR/p6_option1_e2.json"
