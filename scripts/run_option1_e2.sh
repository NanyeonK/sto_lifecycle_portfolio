#!/usr/bin/env bash
# Run E2_2L Option 1 baseline (6D state with per-period tau_buy)
# Key test: does mean_xB > 0 emerge at ell=A? (hedge channel test H1)
# Compare CEV against v3 E1_2L_full (Option 3 baseline = +4.255%)
# Usage:  bash scripts/run_option1_e2.sh
# Server: server1 with Julia at ~/.local/bin/julia

set -e

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

export REGIME=E2_2L
export N_W=15
export N_Z=5
export N_X_PREV=3
export X_PREV_MAX=1.0
export ASSET_GRID_SIZE=9
export GH_NODES=3
export TAU_SELL=0.06
export TAU_BUY=0.025
export TAU_TOKEN=0.01
export P_RELOCATE_WORKING=0.06
export P_RELOCATE_RETIRED=0.02
export RHO_AB=0.50
export GAMMA=5.0
export BETA=0.96
export RF=1.02
export EQUITY_PREMIUM=0.04
export LTV_MAX=0.0
export SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json"

echo "=== v4 E2_2L Option 1 baseline (tau_buy per-period) ==="
echo "State: (t, w, z, ell, x_A_prev, x_B_prev)"
echo "x_prev grid: N=$N_X_PREV, max=$X_PREV_MAX → {0.0, 0.5, 1.0}"
echo "Hypothesis H1: mean_xB > 0 at ell=A (hedge channel activates)"
echo "Started: $(date)"

JULIA_NUM_THREADS=1 \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo "Done: $(date)"
echo "JSON: $OUTDIR/p6_option1_e2.json"
echo ""
echo "Hypotheses to check after run:"
echo "  H1: mean_xB_t1_xprev0_ellA > 0.0  (hedge channel activates)"
echo "  H2: V_t1_midpoint_ellA_xprev0 > v3 E1_2L_full V (-1408.66) → CEV > 4.255%"
echo "  H3: hedge channel CEV ~ 0.5-1.5%"
