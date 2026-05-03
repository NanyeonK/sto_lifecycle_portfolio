#!/usr/bin/env bash
# Run E1_2L baseline with v4 solver (Path B Option 1)
# Execute on server1 from repo root: bash scripts/run_option1_e1.sh
set -euo pipefail

OUT=output/diagnostics
mkdir -p "$OUT"

echo "[run_option1_e1] Starting E1_2L v4 baseline..."
REGIME=E1_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=9 X_GRID_SIZE=5 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 SIGMA_DIV=0.10 \
SUMMARY_JSON_PATH="$OUT/p6_option1_e1.json" \
JULIA_NUM_THREADS=1 \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUT/p6_option1_e1_stdout.log"

echo "[run_option1_e1] Done. JSON: $OUT/p6_option1_e1.json"
