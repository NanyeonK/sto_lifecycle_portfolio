#!/usr/bin/env bash
# v4 E1_2L baseline run.
# State: (t, w, z, ell, x_A_prev, x_B_prev) with N_X_PREV=3.
# Grids: N_W=15, N_Z=5 (reduced vs v3 to offset 9x state expansion).
# Run on server1: bash scripts/run_option1_e1.sh
set -euo pipefail

mkdir -p output/diagnostics

REGIME=E1_2L \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
N_W=15 \
N_Z=5 \
ASSET_GRID_SIZE=7 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json \
julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e1_stdout.log

echo "E1_2L v4 complete."
