#!/usr/bin/env bash
# Run E2_2L_V4 baseline (v4 Option 1 state extension)
# Usage: bash scripts/run_option1_e2.sh
# Env overrides: N_W, N_Z, N_X_PREV, X_PREV_MAX, ASSET_GRID_SIZE, X_GRID_SIZE

set -euo pipefail
cd "$(dirname "$0")/.."

OUT="output/diagnostics/p6_option1_e2.json"
LOG="output/diagnostics/p6_option1_e2_stdout.log"

echo "=== Option 1 E2_2L_V4 baseline ===" | tee "$LOG"
date | tee -a "$LOG"

REGIME=E2_2L_V4 \
N_W=${N_W:-15} \
N_Z=${N_Z:-5} \
N_X_PREV=${N_X_PREV:-3} \
X_PREV_MAX=${X_PREV_MAX:-1.5} \
ASSET_GRID_SIZE=${ASSET_GRID_SIZE:-7} \
X_GRID_SIZE=${X_GRID_SIZE:-4} \
GH_NODES=${GH_NODES:-3} \
SUMMARY_JSON_PATH="$OUT" \
julia src/vfi_solver_v4.jl 2>&1 | tee -a "$LOG"

echo "Done. Summary: $OUT" | tee -a "$LOG"
