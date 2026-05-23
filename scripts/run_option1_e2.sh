#!/usr/bin/env bash
# Run E2_2L baseline with v4 solver (Option 1 full state extension).
# Same grid settings as run_option1_e1.sh for direct CEV comparison.
# Run on server1 inside tmux session sto_lifecycle_portfolio.

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

REGIME=E2_2L
SUMMARY="$OUTDIR/p6_option1_e2.json"
LOG="$OUTDIR/p6_option1_e2_stdout.log"

echo "=== Option 1 E2_2L baseline ===" | tee "$LOG"
date | tee -a "$LOG"

REGIME=$REGIME \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
SUMMARY_JSON_PATH="$SUMMARY" \
julia src/vfi_solver_v4.jl 2>&1 | tee -a "$LOG"

echo "Done. Summary: $SUMMARY" | tee -a "$LOG"
