#!/usr/bin/env bash
# run_option1_e1.sh — E1_2L baseline with v4 6D state (Option 1)
# Run on server1:  bash scripts/run_option1_e1.sh
# Expected wall time: ~2-3 hours single thread

set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$REPO/output/diagnostics"

export REGIME=E1_2L
export N_W=15
export N_Z=5
export N_X_PREV=3
export X_PREV_MAX=1.0
export ASSET_GRID_SIZE=9
export GH_NODES=3
export TAU_SELL=0.06
export TAU_BUY=0.025
export TAU_TOKEN=0.005
export P_RELOCATE_WORKING=0.06
export P_RELOCATE_RETIRED=0.02
export RHO_AB=0.50
export GAMMA=5.0
export BETA=0.96
export RHO=0.05
export M=0.01
export SIGMA_DIV=0.10
export SUMMARY_JSON_PATH="$REPO/output/diagnostics/p6_option1_e1.json"

echo "=== v4 E1_2L Option 1 ==="
date
JULIA_NUM_THREADS=1 julia "$REPO/src/vfi_solver_v4.jl" 2>&1 | tee "$REPO/output/diagnostics/p6_option1_e1_stdout.log"
echo "Done: $(date)"
