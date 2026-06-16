#!/usr/bin/env bash
# Run E2_2L under v4 Option 1 full state extension.
# Continuous fractional tokens with proper tau_buy on positive deltas.
# Pre-holding x_B at ell=A saves tau_buy on relocation to B.
# Outputs: output/diagnostics/p6_option1_e2.json + stdout log.
# Usage: bash scripts/run_option1_e2.sh

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"
JSON_PATH="$OUTDIR/p6_option1_e2.json"
LOG_PATH="$OUTDIR/p6_option1_e2_stdout.log"

echo "[run_option1_e2] starting E2_2L v4 (Option 1 full state) at $(date)"

REGIME=E2_2L \
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

echo "[run_option1_e2] done at $(date). JSON: $JSON_PATH"
