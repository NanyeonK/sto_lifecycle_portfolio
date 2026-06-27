#!/usr/bin/env bash
# Run E1_2L baseline with v4 solver (Option 1: proper x_prev state extension).
# N_X_PREV=3 ({0.0,0.75,1.5}), N_W=15, N_Z=5.
# Note: X_PREV_MAX=1.5 so grid is {0.0, 0.75, 1.5}. x=1.0 is NOT exactly on this grid.
# For E1_2L we need x=1.0 on the grid. Use N_X_PREV=3, X_PREV_MAX=1.0 → {0.0, 0.5, 1.0}.
set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p output/diagnostics

REGIME=E1_2L \
N_W=15 \
N_Z=5 \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 \
GH_NODES=3 \
GAMMA=5.0 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.01 \
RHO_AB=0.50 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json \
  julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e1_stdout.log

echo "E1_2L v4 run complete."
