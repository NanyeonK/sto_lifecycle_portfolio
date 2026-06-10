#!/usr/bin/env bash
# run_option1_e1.sh — E1_2L baseline under v4 6D state extension (Option 1)
#
# Usage: bash scripts/run_option1_e1.sh
# Output: output/diagnostics/p6_option1_e1.json + stdout log
#
# Grid: N_W=15, N_Z=5, N_X_PREV=3 (x_prev ∈ {0.0, 0.5, 1.0}), ASSET=9
#   → ~4.6x v3 compute; expect ~2-3h wall on server1 (single thread).
#
# X_PREV_MAX=1.0 so grid endpoint maps exactly to x=1 (binary own in E1_2L).
# For E1_2L only grid[1]=0.0 and grid[end]=1.0 are used as choices.

set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p output/diagnostics

export REGIME=E1_2L
export N_W=15
export N_Z=5
export N_X_PREV=3
export X_PREV_MAX=1.0
export ASSET_GRID_SIZE=9
export GH_NODES=3
export GAMMA=5.0
export BETA=0.96
export RF=1.02
export EQUITY_PREMIUM=0.04
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
export SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json

echo "[$(date)] Starting E1_2L v4 baseline (Option 1)..."
julia --project=. src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e1_stdout.log
echo "[$(date)] Done. Results → output/diagnostics/p6_option1_e1.json"
