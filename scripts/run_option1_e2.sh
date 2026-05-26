#!/usr/bin/env bash
# run_option1_e2.sh — Run E2_2L baseline with v4 (Option 1: proper tau_buy on deltas)
# Usage: bash scripts/run_option1_e2.sh
# Output: output/diagnostics/p6_option1_e2.json
#
# Key hypothesis checks (from next_actions.md):
#   H1: mean_xB > 0 at ell=A (cross-location hedge activates)
#   H2: CEV(E2_2L_v4 vs E1_2L_v4) > 4.255% (Option 3 baseline)
#   H3: Hedge channel CEV(E2_2L_v4 vs E2_2L_v3) ≈ 0.5-1.5%

set -euo pipefail
mkdir -p output/diagnostics

REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 SIGMA_DIV=0.10 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json \
julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e2_stdout.log

echo "E2_2L v4 run complete. Check output/diagnostics/p6_option1_e2.json"
echo ""
echo "After both runs complete, check H1: mean_xB_t1_xprev00_ellA > 0 in e2.json"
