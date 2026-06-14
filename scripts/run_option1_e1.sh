#!/usr/bin/env bash
# E1_2L baseline at v4 settings (Option 1 state extension)
# Run on server1: bash scripts/run_option1_e1.sh
# Output: output/diagnostics/p6_option1_e1.json

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

export REGIME=E1_2L
export N_X_PREV=3
export X_PREV_MAX=1.0
export N_W=15
export N_Z=5
export ASSET_GRID_SIZE=7
export GH_NODES=3
export RHO=0.05
export M=0.01
export SIGMA_H=0.115
export SIGMA_DIV=0.10
export RHO_AB=0.50
export P_RELOCATE_WORKING=0.06
export P_RELOCATE_RETIRED=0.02
export TAU_SELL=0.06
export TAU_BUY=0.025
export TAU_TOKEN=0.01
export LTV_MAX=0.0
export SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e1.json"

echo "[$(date)] v4 E1_2L baseline starting (Option 1)"
time JULIA_NUM_THREADS=1 julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e1_stdout.log"
echo "[$(date)] v4 E1_2L baseline done. JSON: $SUMMARY_JSON_PATH"
