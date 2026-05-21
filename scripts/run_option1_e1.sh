#!/usr/bin/env bash
# run_option1_e1.sh — E1_2L baseline at v4 grid settings
# Produces comparison baseline for CEV(E2_2L_v4 vs E1_2L)
# Run from repo root: bash scripts/run_option1_e1.sh
set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

JULIA="${JULIA:-julia}"
SOLVER="src/vfi_solver_v4.jl"

echo "=== Option 1: E1_2L baseline (v4 grid, 4D state) ==="
echo "  Grid: N_W=15, N_Z=5, ASSET=7, GH=3"
echo "  tau_sell=0.06, tau_buy=0.025 (apply_at_reloc=1)"
date

JULIA_NUM_THREADS=1 \
REGIME=E1_2L \
N_W=15 N_Z=5 ASSET_GRID_SIZE=7 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 APPLY_TAU_BUY=1 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
N_X_PREV=3 X_PREV_MAX=1.5 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e1.json" \
"$JULIA" "$SOLVER"

echo "Done. Output: $OUTDIR/p6_option1_e1.json"
date
