#!/usr/bin/env bash
# run_option1_e2.sh — E2_2L baseline at v4 6D state settings
# Usage: bash scripts/run_option1_e2.sh
# Runs on server1; expects Julia in PATH and JSON3 available.
set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"
LOG="$OUTDIR/p6_option1_e2_stdout.log"

REGIME=E2_2L \
N_W=15 \
N_Z=5 \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 \
GH_NODES=3 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
JULIA_NUM_THREADS=1 \
julia src/vfi_solver_v4.jl 2>&1 | tee "$LOG"

echo "E2_2L done. Summary at $OUTDIR/p6_option1_e2.json"
