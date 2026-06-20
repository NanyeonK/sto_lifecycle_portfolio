#!/usr/bin/env bash
# run_option1_e1.sh — Run E1_2L baseline with v4 solver (Option 1 state extension)
#
# Usage:
#   bash scripts/run_option1_e1.sh            # default coarse grids (N_W=15, N_Z=5, N_X_PREV=3)
#   FULL=1 bash scripts/run_option1_e1.sh     # full grids (N_W=40, N_Z=9)
#
# Output:
#   output/diagnostics/p6_option1_e1.json     — summary JSON
#   output/diagnostics/p6_option1_e1_stdout.log

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$REPO_ROOT/output/diagnostics"
mkdir -p "$OUT_DIR"

FULL="${FULL:-0}"

if [ "$FULL" = "1" ]; then
    export N_W=40 N_Z=9 N_X_PREV=3 X_PREV_MAX=1.0 ASSET_GRID_SIZE=13 GH_NODES=3
    LABEL="full"
else
    export N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 ASSET_GRID_SIZE=7 GH_NODES=3
    LABEL="coarse"
fi

export REGIME=E1_2L
export TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01
export RHO_AB=0.50 P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02
export GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04
export RHO=0.05 M=0.01 SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016
export LTV_MAX=0.0
export SUMMARY_JSON_PATH="$OUT_DIR/p6_option1_e1_${LABEL}.json"

echo "=== v4 E1_2L Option 1 ($LABEL grids) ==="
echo "  N_W=$N_W N_Z=$N_Z N_X_PREV=$N_X_PREV ASSET_GRID_SIZE=$ASSET_GRID_SIZE"
echo "  tau_sell=$TAU_SELL tau_buy=$TAU_BUY tau_token=$TAU_TOKEN"
echo "  rho_AB=$RHO_AB p_relocate_working=$P_RELOCATE_WORKING"
date

julia --project="$REPO_ROOT" "$REPO_ROOT/src/vfi_solver_v4.jl" \
    2>&1 | tee "$OUT_DIR/p6_option1_e1_${LABEL}_stdout.log"

echo "Done. Summary at: $SUMMARY_JSON_PATH"
