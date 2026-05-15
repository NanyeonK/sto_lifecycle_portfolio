#!/usr/bin/env bash
# E1_2L baseline run for vfi_solver_v4.jl (Option 1 state extension)
# Coarse grids: N_W=15, N_Z=5, N_X_PREV=3, X_PREV_MAX=1.0
# Expected wall time: ~2-3 hours single thread (4.6x v3 baseline estimate)
#
# Usage:
#   bash scripts/run_option1_e1.sh
#   JULIA_NUM_THREADS=4 bash scripts/run_option1_e1.sh   # parallel (experimental)

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$REPO/output/diagnostics"
mkdir -p "$OUT"

LOG="$OUT/p6_option1_e1_stdout.log"
JSON="$OUT/p6_option1_e1.json"

echo "=== Option 1 E1_2L run ==="
echo "Branch: $(git -C "$REPO" rev-parse --abbrev-ref HEAD)"
echo "Commit: $(git -C "$REPO" rev-parse --short HEAD)"
echo "Output: $JSON"
echo ""

REGIME=E1_2L \
N_W=15 \
N_Z=5 \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.005 \
RHO_AB=0.5 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
SUMMARY_JSON_PATH="$JSON" \
julia "$REPO/src/vfi_solver_v4.jl" 2>&1 | tee "$LOG"

echo ""
echo "Summary JSON: $JSON"
echo "Log:          $LOG"
