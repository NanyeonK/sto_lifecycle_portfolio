#!/usr/bin/env bash
# Run E2_2L Option 1 baseline at v4 6D state settings.
set -euo pipefail

mkdir -p output/diagnostics

REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json \
julia --project=. src/vfi_solver_v4.jl \
  2>&1 | tee output/diagnostics/p6_option1_e2_stdout.log

echo "E2_2L v4 run complete."
