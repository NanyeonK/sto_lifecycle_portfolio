#!/usr/bin/env bash
# Run E1_2L baseline with v4 solver (Option 1 full state extension).
# Expected wall time: ~2.5 hours single thread at N_X_PREV=3.
# Output: output/diagnostics/p6_option1_e1.json
set -euo pipefail

mkdir -p output/diagnostics

REGIME=E1_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.5 GAMMA=5.0 BETA=0.96 RF=1.02 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json \
julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e1_stdout.log
