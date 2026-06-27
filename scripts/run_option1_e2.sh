#!/usr/bin/env bash
# Run E2_2L baseline with v4 solver (Option 1: proper x_prev state extension).
# x_prev grid: {0.0, 0.5, 1.0} — allows pre-holding B at half/full unit before relocation.
set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p output/diagnostics

REGIME=E2_2L \
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
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json \
  julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e2_stdout.log

echo "E2_2L v4 run complete."
