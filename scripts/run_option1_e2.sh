#!/usr/bin/env bash
# Run E2_2L baseline with v4 (6D state, Option 1 tau_buy on positive deltas).
# Usage: bash scripts/run_option1_e2.sh
# Writes output to output/diagnostics/p6_option1_e2.json

set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p output/diagnostics

REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=9 GH_NODES=3 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
RHO_AB=0.50 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json \
julia src/vfi_solver_v4.jl \
  2>&1 | tee output/diagnostics/p6_option1_e2_stdout.log

echo "E2_2L v4 done. Summary: output/diagnostics/p6_option1_e2.json"
