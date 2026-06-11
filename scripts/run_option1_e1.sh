#!/usr/bin/env bash
# E1_2L baseline run for v4 (Option 1 full state extension).
# N_W=15, N_Z=5, N_X_PREV=3 — coarse grid, ~2-3h wall single thread.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"
mkdir -p output/diagnostics

echo "=== v4 E1_2L baseline ==="
REGIME=E1_2L \
N_W=15 N_Z=5 \
N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
APPLY_TAU_BUY=0 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json \
JULIA_NUM_THREADS=1 \
julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e1_stdout.log

echo "Done. Summary: output/diagnostics/p6_option1_e1.json"
