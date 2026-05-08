#!/usr/bin/env bash
# Run E2_2L baseline at v4 Option 1 settings (6D state, proper tau_buy).
# Expected wall time: ~2-3 hours single thread at N_W=15, N_Z=5, N_X_PREV=3.
# Run from repo root: bash scripts/run_option1_e2.sh
set -euo pipefail

OUT=output/diagnostics
mkdir -p "$OUT"

JULIA_NUM_THREADS=1 \
  REGIME=E2_2L \
  N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
  ASSET_GRID_SIZE=7 GH_NODES=3 \
  TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
  P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
  RHO_AB=0.50 \
  SUMMARY_JSON_PATH="$OUT/p6_option1_e2.json" \
  julia src/vfi_solver_v4.jl 2>&1 | tee "$OUT/p6_option1_e2_stdout.log"

echo "E2_2L v4 done. Summary: $OUT/p6_option1_e2.json"
