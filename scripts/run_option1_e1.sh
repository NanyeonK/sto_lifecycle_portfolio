#!/usr/bin/env bash
# Run E1_2L baseline with v4 solver (Option 1 full state extension).
# Usage:  bash scripts/run_option1_e1.sh
# Output: output/diagnostics/p6_option1_e1.json  +  stdout log

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

LOG="$OUTDIR/p6_option1_e1_stdout.log"
JSON="$OUTDIR/p6_option1_e1.json"

echo "=== v4 E1_2L Option 1 baseline ==="
echo "Log: $LOG"
echo "JSON: $JSON"

REGIME=E1_2L \
SUMMARY_JSON_PATH="$JSON" \
N_W=15 \
N_Z=5 \
N_X_PREV=3 \
X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=7 \
X_GRID_SIZE=4 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
LTV_MAX=0.0 \
JULIA_NUM_THREADS=1 \
  julia src/vfi_solver_v4.jl 2>&1 | tee "$LOG"

echo "=== E1_2L done. JSON written to $JSON ==="
