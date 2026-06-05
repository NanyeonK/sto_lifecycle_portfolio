#!/usr/bin/env bash
# Run E2_2L baseline with v4 solver (Option 1 full state extension)
# Coarse grid: ~2-3 h wall on server1 single thread.
# Usage: bash scripts/run_option1_e2.sh
set -euo pipefail
mkdir -p output/diagnostics
REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=5 GH_NODES=3 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json \
julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e2_stdout.log
echo "E2_2L v4 done — results in output/diagnostics/p6_option1_e2.json"
