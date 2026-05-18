#!/usr/bin/env bash
# run_option1_e2.sh — E2_2L baseline for v4 Option 1 state extension
# Usage: bash scripts/run_option1_e2.sh
# Output: output/diagnostics/p6_option1_e2.json + stdout log

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
mkdir -p "$REPO_ROOT/output/diagnostics"

echo "[run_option1_e2] Starting E2_2L v4 at $(date)"

JULIA_NUM_THREADS=1 \
REGIME=E2_2L \
N_W=15 \
N_Z=5 \
N_X_PREV=3 \
X_PREV_MAX=2.0 \
ASSET_GRID_SIZE=7 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
SUMMARY_JSON_PATH="$REPO_ROOT/output/diagnostics/p6_option1_e2.json" \
julia "$REPO_ROOT/src/vfi_solver_v4.jl" \
    2>&1 | tee "$REPO_ROOT/output/diagnostics/p6_option1_e2_stdout.log"

echo "[run_option1_e2] Done at $(date)"
echo "[run_option1_e2] Summary JSON: $REPO_ROOT/output/diagnostics/p6_option1_e2.json"
