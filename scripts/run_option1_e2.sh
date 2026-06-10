#!/usr/bin/env bash
# run_option1_e2.sh — E2_2L baseline under v4 6D state extension (Option 1)
#
# Usage: bash scripts/run_option1_e2.sh
# Output: output/diagnostics/p6_option1_e2.json + stdout log
#
# Grid: N_W=15, N_Z=5, N_X_PREV=3, ASSET=9.
# X_PREV_MAX=1.0 → x_prev ∈ {0.0, 0.5, 1.0} for both locations.
# E2_2L uses all 9 (x_A_new, x_B_new) ∈ grid^2 combinations.
#
# Key test: does mean_xB_t1_entry_ellA > 0?
#   YES → hedge channel activates (pre-holding B while at A)
#   NO  → mechanism still dormant; escalate to Option C or accept Path D
#
# After both runs complete, compute CEV:
#   CEV = [(V_e2/V_e1)^(1/(1-gamma)) - 1] at entry state (ix_A=1, ix_B=1)

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
export SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json

echo "[$(date)] Starting E2_2L v4 baseline (Option 1)..."
julia --project=. src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e2_stdout.log
echo "[$(date)] Done. Results → output/diagnostics/p6_option1_e2.json"
