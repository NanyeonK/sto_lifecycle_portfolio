#!/usr/bin/env bash
# Run E2_2L baseline with v4 (Option 1 state extension).
# Usage: bash scripts/run_option1_e2.sh
# Expected wall: ~2-3 h on server1 single thread.
# Run AFTER run_option1_e1.sh so CEV comparison is ready.

set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p output/diagnostics

export REGIME=E2_2L
export N_W=15
export N_Z=5
export N_X_PREV=3
export X_PREV_MAX=1.0
export ASSET_GRID_SIZE=9
export GH_NODES=3
export TAU_SELL=0.06
export TAU_BUY=0.025
export TAU_TOKEN=0.01
export RHO_AB=0.50
export P_RELOCATE_WORKING=0.06
export P_RELOCATE_RETIRED=0.02
export SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json

JULIA_NUM_THREADS=1 julia src/vfi_solver_v4.jl \
    2>&1 | tee output/diagnostics/p6_option1_e2_stdout.log

echo "E2_2L v4 done. JSON at output/diagnostics/p6_option1_e2.json"
