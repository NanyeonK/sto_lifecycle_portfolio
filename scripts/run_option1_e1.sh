#!/usr/bin/env bash
# Run E1_2L baseline under v4 (Option 1 full state extension).
# Coarse grids: N_W=15, N_Z=5, N_X_PREV=3, X_PREV_MAX=1.0 (binary {0,0.5,1}).
# Expected wall time: ~2-3 hours on server1 (single thread).
# Usage: bash scripts/run_option1_e1.sh

set -euo pipefail

OUT_DIR="output/diagnostics"
mkdir -p "$OUT_DIR"

JULIA="${JULIA:-julia}"
SOLVER="src/vfi_solver_v4.jl"
OUTFILE="$OUT_DIR/p6_option1_e1.json"
LOGFILE="$OUT_DIR/p6_option1_e1_stdout.log"

echo "[$(date)] Running E1_2L v4 baseline..." | tee "$LOGFILE"

REGIME=E1_2L \
N_W=15 N_Z=5 \
N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 \
GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5.0 BETA=0.96 \
SUMMARY_JSON_PATH="$OUTFILE" \
JULIA_NUM_THREADS=1 \
  "$JULIA" "$SOLVER" 2>&1 | tee -a "$LOGFILE"

echo "[$(date)] E1_2L v4 DONE — output: $OUTFILE" | tee -a "$LOGFILE"
