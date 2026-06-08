#!/usr/bin/env bash
# run_option1_e1.sh — E1_2L baseline run for v4 Option 1 solver
# State: (t, w, z, ell, x_A_prev, x_B_prev) — 6D
# tau_buy applied per-period on positive delta at current location.
#
# Reduced grids to compensate for 9x state-space factor:
#   N_W=15, N_Z=5, N_X_PREV=3 (X_PREV_MAX=2.0 → {0.0, 1.0, 2.0})
# Expected wall time: ~2–2.5 h single thread.
#
# Usage:  bash scripts/run_option1_e1.sh [output_dir]
# Output: output/diagnostics/p6_option1_e1.json  (+ stdout log)

set -euo pipefail

OUTPUT_DIR="${1:-output/diagnostics}"
mkdir -p "$OUTPUT_DIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="$OUTPUT_DIR/p6_option1_e1_stdout_${TIMESTAMP}.log"
JSON_FILE="$OUTPUT_DIR/p6_option1_e1.json"

echo "=== v4 E1_2L Option 1 baseline run ===" | tee "$LOG_FILE"
echo "  start: $(date)"                          | tee -a "$LOG_FILE"
echo "  log:   $LOG_FILE"                        | tee -a "$LOG_FILE"
echo "  json:  $JSON_FILE"                       | tee -a "$LOG_FILE"
echo ""                                           | tee -a "$LOG_FILE"

JULIA_NUM_THREADS=1 \
  REGIME=E1_2L \
  N_W=15 N_Z=5 \
  N_X_PREV=3 X_PREV_MAX=2.0 \
  ASSET_GRID_SIZE=9 X_GRID_SIZE=5 GH_NODES=3 \
  TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
  P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
  RHO_AB=0.50 \
  LTV_MAX=0.0 \
  SUMMARY_JSON_PATH="$JSON_FILE" \
  julia src/vfi_solver_v4.jl 2>&1 | tee -a "$LOG_FILE"

echo "" | tee -a "$LOG_FILE"
echo "=== done: $(date) ===" | tee -a "$LOG_FILE"
