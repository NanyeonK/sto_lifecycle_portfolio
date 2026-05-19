#!/usr/bin/env bash
# run_option1_e2.sh — E2_2L Option 1 run (6D state, proper tau_buy hedge)
# Run on server1. Output: output/diagnostics/p6_option1_e2.json
# Usage: bash scripts/run_option1_e2.sh

set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p output/diagnostics

# Key hypothesis: mean_xB > 0 at ell=A (hedge mechanism activates).
# x_prev grid {0, 0.5, 1.0}: pre-holding x_B costs tau_buy; keeping across move costs 0.
REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5.0 BETA=0.96 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json \
julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e2_stdout.log

echo "E2_2L v4 done. JSON: output/diagnostics/p6_option1_e2.json"
echo "Check mean_xB_t1_ellA > 0 to verify hedge mechanism activates."
