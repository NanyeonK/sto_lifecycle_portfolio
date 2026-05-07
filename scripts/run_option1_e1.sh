#!/usr/bin/env bash
# run_option1_e1.sh — E1_2L baseline run for Option 1 (6D state, v4 solver)
# Run on server1 from repo root: bash scripts/run_option1_e1.sh
#
# Expected wall time: ~2-3h single thread (N_W=15, N_Z=5, N_X_PREV=3, GH=3)
# Output: output/diagnostics/p6_option1_e1.json  + stdout log

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="$REPO_ROOT/output/diagnostics"
mkdir -p "$OUT_DIR"

LOG="$OUT_DIR/p6_option1_e1_stdout.log"
JSON="$OUT_DIR/p6_option1_e1.json"

echo "=== Option 1 E1_2L baseline (v4 solver) ===" | tee "$LOG"
date | tee -a "$LOG"

REGIME=E1_2L \
N_W=15 \
N_Z=5 \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5.0 \
BETA=0.96 \
SUMMARY_JSON_PATH="$JSON" \
  julia "$REPO_ROOT/src/vfi_solver_v4.jl" 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== E1_2L run complete ===" | tee -a "$LOG"
date | tee -a "$LOG"
echo "Summary JSON: $JSON" | tee -a "$LOG"
