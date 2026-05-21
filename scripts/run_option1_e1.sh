#!/usr/bin/env bash
# run_option1_e1.sh — E1_2L baseline run with v4 Option 1 state extension
# Usage: bash scripts/run_option1_e1.sh
# Prereq: Julia installed, JSON3 available in global Julia env.
#
# Grid: N_W=15, N_Z=5, N_X_PREV=3, GH_NODES=3 (smoke-scale default)
# Expected runtime: ~2-3 hours per regime on server1 (single thread).

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "[$(date)] Starting E1_2L v4 Option 1 run"

REGIME=E1_2L \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
N_W=15 \
N_Z=5 \
ASSET_GRID_SIZE=9 \
GH_NODES=3 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.005 \
APPLY_TAU_BUY=1 \
RHO_AB=0.50 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e1.json" \
  julia src/vfi_solver_v4.jl \
  2>&1 | tee "$OUTDIR/p6_option1_e1_stdout.log"

echo "[$(date)] E1_2L v4 done. Output: $OUTDIR/p6_option1_e1.json"
