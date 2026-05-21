#!/usr/bin/env bash
# run_option1_e2.sh — E2_2L_v4 (Option 1: 6D state, delta tx_cost)
# Key hypothesis: mean_xB > 0 at ell=A (hedge channel activates).
# Run from repo root: bash scripts/run_option1_e2.sh
set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

JULIA="${JULIA:-julia}"
SOLVER="src/vfi_solver_v4.jl"

echo "=== Option 1: E2_2L_v4 (6D state, delta tx_cost) ==="
echo "  Grid: N_W=15, N_Z=5, ASSET=7, N_X_PREV=3"
echo "  tau_buy=0.025 (per-period delta), tau_token=0.005"
echo "  Expected: mean_xB > 0 at ell=A (pre-buy hedge active)"
date

JULIA_NUM_THREADS=1 \
REGIME=E2_2L \
N_W=15 N_Z=5 ASSET_GRID_SIZE=7 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
N_X_PREV=3 X_PREV_MAX=1.5 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
"$JULIA" "$SOLVER"

echo "Done. Output: $OUTDIR/p6_option1_e2.json"
date
