#!/bin/bash
# run_option1_e1.sh — E1_2L baseline for v4 Option 1 full state extension
# Run on server1; ~2-3 hours wall single thread (N_W=15, N_Z=5, N_X_PREV=3)
#
# Usage:
#   bash scripts/run_option1_e1.sh
#   JULIA_NUM_THREADS=4 bash scripts/run_option1_e1.sh  # multi-thread

set -euo pipefail

OUTPUT_DIR="output/diagnostics"
mkdir -p "$OUTPUT_DIR"

SUMMARY_JSON="$OUTPUT_DIR/p6_option1_e1.json"
STDOUT_LOG="$OUTPUT_DIR/p6_option1_e1_stdout.log"

echo "=== Option 1 E1_2L run: $(date) ==="
echo "Output: $SUMMARY_JSON"
echo "Log:    $STDOUT_LOG"

REGIME=E1_2L \
N_W=15 \
N_Z=5 \
N_X_PREV=3 \
X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=7 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5.0 \
BETA=0.96 \
SUMMARY_JSON_PATH="$SUMMARY_JSON" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$STDOUT_LOG"

echo "=== E1_2L done: $(date) ==="
echo "Summary at: $SUMMARY_JSON"
