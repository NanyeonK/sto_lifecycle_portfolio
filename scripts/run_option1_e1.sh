#!/usr/bin/env bash
# Run E1_2L baseline under v4 Option 1 state extension.
# Binary ownership, tau_buy charged on buying x_ell (not via approximation).
# Outputs: output/diagnostics/p6_option1_e1.json + stdout log.
# Usage: bash scripts/run_option1_e1.sh

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"
JSON_PATH="$OUTDIR/p6_option1_e1.json"
LOG_PATH="$OUTDIR/p6_option1_e1_stdout.log"

echo "[run_option1_e1] starting E1_2L v4 at $(date)"

REGIME=E1_2L \
N_X_PREV=3 \
X_PREV_MAX=2.0 \
N_W=15 \
N_Z=5 \
ASSET_GRID_SIZE=9 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
SIGMA_DIV=0.10 \
GAMMA=5.0 \
BETA=0.96 \
SUMMARY_JSON_PATH="$JSON_PATH" \
    julia src/vfi_solver_v4.jl 2>&1 | tee "$LOG_PATH"

echo "[run_option1_e1] done at $(date). JSON: $JSON_PATH"
