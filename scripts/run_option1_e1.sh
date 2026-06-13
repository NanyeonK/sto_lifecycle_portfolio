#!/usr/bin/env bash
# Run E1_2L baseline for Option 1 (v4 solver, 6D state)
# Usage: bash scripts/run_option1_e1.sh
# Output: output/diagnostics/p6_option1_e1.json

set -euo pipefail
cd "$(dirname "$0")/.."

REGIME=E1_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json \
  julia src/vfi_solver_v4.jl \
  | tee output/diagnostics/p6_option1_e1_stdout.log

echo "E1_2L done. JSON: output/diagnostics/p6_option1_e1.json"
