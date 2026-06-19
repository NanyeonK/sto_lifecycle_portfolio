#!/usr/bin/env bash
# Run E2_2L baseline under v4 (Option 1 full state extension).
# x_prev_grid = {0.0, 0.5, 1.0} — 3 points covering renter/partial/full-owner.
# Hypothesis: mean_xB > 0 at ell=A will emerge with proper delta tx_cost.
#
# Usage (on server1):
#   bash scripts/run_option1_e2.sh
#   JULIA_NUM_THREADS=4 N_X_PREV=4 X_PREV_MAX=1.5 bash scripts/run_option1_e2.sh

set -e

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

export REGIME="E2_2L"
export N_W="${N_W:-15}"
export N_Z="${N_Z:-5}"
export N_X_PREV="${N_X_PREV:-3}"
export X_PREV_MAX="${X_PREV_MAX:-1.0}"
export ASSET_GRID_SIZE="${ASSET_GRID_SIZE:-9}"
export GH_NODES="${GH_NODES:-3}"
export TAU_SELL="${TAU_SELL:-0.06}"
export TAU_BUY="${TAU_BUY:-0.025}"
export TAU_TOKEN="${TAU_TOKEN:-0.01}"
export LTV_MAX="${LTV_MAX:-0.0}"
export SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json"

echo "=== v4 E2_2L baseline (Option 1) ==="
echo "    N_W=$N_W  N_Z=$N_Z  N_X_PREV=$N_X_PREV  X_PREV_MAX=$X_PREV_MAX"
echo "    tau_sell=$TAU_SELL  tau_buy=$TAU_BUY  tau_token=$TAU_TOKEN"
echo "    Output: $OUTDIR/p6_option1_e2.json"
echo "    Key check: mean_xB_t1_feasible_ellA > 0 (hedge motive activates)"
echo ""

time julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo ""
echo "Done. Summary at $OUTDIR/p6_option1_e2.json"
echo "Key diagnostic: grep mean_xB $OUTDIR/p6_option1_e2.json"
