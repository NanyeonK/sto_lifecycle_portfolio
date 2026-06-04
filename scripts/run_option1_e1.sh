#!/usr/bin/env bash
# Run E1_2L baseline with Option 1 (v4 solver, full state extension)
# Usage: bash scripts/run_option1_e1.sh
set -euo pipefail

REGIME=E1_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
RHO_AB=0.50 P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
GAMMA=5.0 BETA=0.96 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json \
  julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e1_stdout.log

echo "E1_2L Option 1 done. JSON: output/diagnostics/p6_option1_e1.json"
