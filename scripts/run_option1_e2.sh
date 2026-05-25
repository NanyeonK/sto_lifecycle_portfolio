#!/usr/bin/env bash
# Run E2_2L baseline with v4 solver (6D state, Option 1 tau_buy on deltas)
# Usage: bash scripts/run_option1_e2.sh
set -e
mkdir -p output/diagnostics

REGIME=E2_2L \
N_X_PREV=3 \
X_PREV_MAX=1.5 \
N_W=15 \
N_Z=5 \
ASSET_GRID_SIZE=7 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.5 \
GAMMA=5.0 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json \
julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e2_stdout.log

echo "E2_2L v4 done. Summary at output/diagnostics/p6_option1_e2.json"
