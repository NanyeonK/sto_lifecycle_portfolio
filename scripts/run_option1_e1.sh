#!/usr/bin/env bash
# run_option1_e1.sh — E1_2L baseline run for v4 Option 1 state extension
# Usage: bash scripts/run_option1_e1.sh
#
# Reduced grids: N_W=15, N_Z=5, N_X_PREV=3 to keep runtime ~2-3h on server1 (single thread).
# Set JULIA_NUM_THREADS=1 explicitly to prevent accidental parallelism.
# Output: output/diagnostics/p6_option1_e1.json + stdout log.

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "[run_option1_e1] Starting E1_2L baseline (v4 Option 1)  $(date)"

JULIA_NUM_THREADS=1 \
REGIME=E1_2L \
N_W=15 \
N_Z=5 \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5.0 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="${OUTDIR}/p6_option1_e1.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "${OUTDIR}/p6_option1_e1_stdout.log"

echo "[run_option1_e1] Done  $(date)"
echo "JSON output: ${OUTDIR}/p6_option1_e1.json"
