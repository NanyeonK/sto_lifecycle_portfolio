#!/usr/bin/env bash
# Run E1_2L baseline at v4 settings (6D state, coarse x_prev grid)
# Compare V values against v3 E1_2L_full to verify structural equivalence.
# Usage:  bash scripts/run_option1_e1.sh
# Server: server1 with Julia at ~/.local/bin/julia

set -e

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

export REGIME=E1_2L
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
export SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e1.json"

echo "=== v4 E1_2L baseline (Option 1) ==="
echo "State: (t, w, z, ell, x_A_prev, x_B_prev)"
echo "x_prev grid: N=$N_X_PREV, max=$X_PREV_MAX"
echo "Started: $(date)"

JULIA_NUM_THREADS=1 \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e1_stdout.log"

echo "Done: $(date)"
echo "JSON: $OUTDIR/p6_option1_e1.json"
