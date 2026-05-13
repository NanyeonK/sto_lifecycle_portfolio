#!/usr/bin/env bash
# Run E2_2L Option 1 baseline (v4 solver, 6D state)
# Default: N_W=15, N_Z=5, N_X_PREV=3, ASSET_GRID_SIZE=5 (~2.5h single thread)
# Override: N_W=5 N_Z=3 ASSET_GRID_SIZE=3 for quick smoke (~5-10 min)
set -e

mkdir -p output/diagnostics

REGIME=E2_2L \
N_W=${N_W:-15} \
N_Z=${N_Z:-5} \
N_X_PREV=${N_X_PREV:-3} \
X_PREV_MAX=${X_PREV_MAX:-1.0} \
ASSET_GRID_SIZE=${ASSET_GRID_SIZE:-5} \
GH_NODES=${GH_NODES:-3} \
TAU_SELL=${TAU_SELL:-0.06} \
TAU_BUY=${TAU_BUY:-0.025} \
TAU_TOKEN=${TAU_TOKEN:-0.005} \
P_RELOCATE_WORKING=${P_RELOCATE_WORKING:-0.06} \
P_RELOCATE_RETIRED=${P_RELOCATE_RETIRED:-0.02} \
RHO_AB=${RHO_AB:-0.50} \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json \
julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e2_stdout.log

echo "E2_2L v4 done. Summary: output/diagnostics/p6_option1_e2.json"
