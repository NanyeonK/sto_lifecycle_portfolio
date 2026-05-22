#!/bin/bash
# run_option1_e2.sh — E2_2L baseline for v4 Option 1 full state extension
# Run on server1 AFTER E1_2L completes; compare CEV(E2_2L_v4 vs E1_2L_v4)
#
# Usage:
#   bash scripts/run_option1_e2.sh
#   JULIA_NUM_THREADS=4 bash scripts/run_option1_e2.sh

set -euo pipefail

OUTPUT_DIR="output/diagnostics"
mkdir -p "$OUTPUT_DIR"

SUMMARY_JSON="$OUTPUT_DIR/p6_option1_e2.json"
STDOUT_LOG="$OUTPUT_DIR/p6_option1_e2_stdout.log"

echo "=== Option 1 E2_2L run: $(date) ==="
echo "Output: $SUMMARY_JSON"
echo "Log:    $STDOUT_LOG"

REGIME=E2_2L \
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

echo "=== E2_2L done: $(date) ==="
echo "Summary at: $SUMMARY_JSON"
echo ""
echo "Next: compute CEV(E2_2L_v4 vs E1_2L_v4) to test H1+H2+H3 hypotheses."
echo "  H1: mean_xB > 0 at ellA  (hedge mechanism activates)"
echo "  H2: CEV(E2_2L_v4 vs E1_2L_v4) > 4.255%  (Option 3 baseline)"
echo "  H3: CEV(E2_2L_v4 vs E2_2L_v3) ≈ 0.5-1.5%  (incremental hedge channel)"
