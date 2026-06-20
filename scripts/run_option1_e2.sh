#!/usr/bin/env bash
# run_option1_e2.sh — Run E2_2L baseline with v4 solver (Option 1 state extension)
#
# Usage:
#   bash scripts/run_option1_e2.sh            # default coarse grids
#   FULL=1 bash scripts/run_option1_e2.sh     # full grids
#   N_X_PREV=5 bash scripts/run_option1_e2.sh # finer x_prev grid
#
# Output:
#   output/diagnostics/p6_option1_e2.json     — summary JSON
#   output/diagnostics/p6_option1_e2_stdout.log
#
# Key hypothesis to check after run:
#   H1: mean_xB > 0 at ell=A  (hedge mechanism activates)
#   H2: CEV(E2_2L_v4 vs E1_2L_v4) > 4.255%  (Option 3 baseline)
#   H3: hedge channel = CEV(E2_2L_v4 vs E2_2L_v3) ≈ 0.5-1.5%

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$REPO_ROOT/output/diagnostics"
mkdir -p "$OUT_DIR"

FULL="${FULL:-0}"

if [ "$FULL" = "1" ]; then
    export N_W=40 N_Z=9 ASSET_GRID_SIZE=13 GH_NODES=3
    LABEL="full"
else
    export N_W=15 N_Z=5 ASSET_GRID_SIZE=7 GH_NODES=3
    LABEL="coarse"
fi

export N_X_PREV="${N_X_PREV:-3}"
export X_PREV_MAX="${X_PREV_MAX:-1.0}"

export REGIME=E2_2L
export TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01
export RHO_AB=0.50 P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02
export GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04
export RHO=0.05 M=0.01 SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016
export LTV_MAX=0.0
export SUMMARY_JSON_PATH="$OUT_DIR/p6_option1_e2_${LABEL}.json"

echo "=== v4 E2_2L Option 1 ($LABEL grids, N_X_PREV=$N_X_PREV) ==="
echo "  N_W=$N_W N_Z=$N_Z N_X_PREV=$N_X_PREV ASSET_GRID_SIZE=$ASSET_GRID_SIZE"
echo "  tau_sell=$TAU_SELL tau_buy=$TAU_BUY tau_token=$TAU_TOKEN"
echo "  rho_AB=$RHO_AB p_relocate_working=$P_RELOCATE_WORKING"
echo "  [Hedge premium per unit x_B: p_relocate × tau_buy = $(echo "$P_RELOCATE_WORKING * $TAU_BUY" | bc -l)]"
date

julia --project="$REPO_ROOT" "$REPO_ROOT/src/vfi_solver_v4.jl" \
    2>&1 | tee "$OUT_DIR/p6_option1_e2_${LABEL}_stdout.log"

echo "Done. Summary at: $SUMMARY_JSON_PATH"
echo ""
echo "Check H1: mean_xB_t1_ellA > 0  (from JSON)"
echo "Check H2: V_t1_midpoint_ellA > E1_2L value (compute CEV after both runs)"
