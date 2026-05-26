#!/usr/bin/env bash
# run_option1_e1.sh — Run E1_2L baseline with v4 (Option 1 state extension)
# Usage: bash scripts/run_option1_e1.sh
# Output: output/diagnostics/p6_option1_e1.json

set -euo pipefail
mkdir -p output/diagnostics

REGIME=E1_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 SIGMA_DIV=0.10 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json \
julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e1_stdout.log

echo "E1_2L v4 run complete. Check output/diagnostics/p6_option1_e1.json"
