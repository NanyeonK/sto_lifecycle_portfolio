#!/usr/bin/env bash
# Run E1_2L baseline under Option 1 (v4 solver, 6D state extension)
# Usage: bash scripts/run_option1_e1.sh
# Output: output/diagnostics/p6_option1_e1.json + p6_option1_e1_stdout.log
set -euo pipefail

mkdir -p output/diagnostics

export REGIME=E1_2L
export N_W=15
export N_Z=5
export N_X_PREV=3
export X_PREV_MAX=1.5
export ASSET_GRID_SIZE=7
export X_GRID_SIZE=4
export GH_NODES=3
export TAU_SELL=0.06
export TAU_BUY=0.025
export TAU_TOKEN=0.005
export P_RELOCATE_WORKING=0.06
export P_RELOCATE_RETIRED=0.02
export RHO_AB=0.50
export LTV_MAX=0.0
export SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json

echo "=== Option 1 E1_2L run ===" | tee output/diagnostics/p6_option1_e1_stdout.log
echo "  N_W=$N_W N_Z=$N_Z N_X_PREV=$N_X_PREV"  | tee -a output/diagnostics/p6_option1_e1_stdout.log
date | tee -a output/diagnostics/p6_option1_e1_stdout.log

JULIA_NUM_THREADS=1 julia src/vfi_solver_v4.jl \
    2>&1 | tee -a output/diagnostics/p6_option1_e1_stdout.log

echo "=== DONE E1_2L ===" | tee -a output/diagnostics/p6_option1_e1_stdout.log
date | tee -a output/diagnostics/p6_option1_e1_stdout.log
