#!/usr/bin/env bash
# Run E2_2L baseline using v4 solver (Option 1 full state extension).
# Expected wall time: ~2-3h on server1 (single thread, small grids).
# Run from repo root: bash scripts/run_option1_e2.sh

set -euo pipefail

OUT="output/diagnostics"
mkdir -p "$OUT"

SUMMARY_JSON="$OUT/p6_option1_e2.json"
LOG="$OUT/p6_option1_e2_stdout.log"

echo "[run_option1_e2] Starting E2_2L v4 baseline at $(date)"

REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
SIGMA_S=0.157 G_H=0.016 SIGMA_H=0.115 SIGMA_DIV=0.10 \
RHO=0.05 M=0.01 \
RHO_AB=0.50 P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$SUMMARY_JSON" \
JULIA_NUM_THREADS=1 \
  julia src/vfi_solver_v4.jl 2>&1 | tee "$LOG"

echo "[run_option1_e2] Done at $(date). Summary: $SUMMARY_JSON"
