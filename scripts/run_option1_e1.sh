#!/usr/bin/env bash
# run_option1_e1.sh — E1_2L baseline run for v4 Option-1 solver
# Usage: bash scripts/run_option1_e1.sh
# Writes: output/diagnostics/p6_option1_e1.json  (summary JSON)
#         output/diagnostics/p6_option1_e1_stdout.log

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="$REPO_ROOT/output/diagnostics"
mkdir -p "$OUT_DIR"

# Coarse grids per Option-1 spec: N_W=15, N_Z=5, N_XPREV=3
# ~2.5h wall on server1 (single thread)
export REGIME=E1_2L
export N_W=15
export N_Z=5
export N_X_PREV=3
export X_PREV_MAX=2.0
export ASSET_GRID_SIZE=7
export X_GRID_SIZE=4
export GH_NODES=3
export SIGMA_DIV=0.10
export RHO_AB=0.50
export P_RELOCATE_WORKING=0.06
export P_RELOCATE_RETIRED=0.02
export TAU_SELL=0.06
export TAU_BUY=0.025
export TAU_TOKEN=0.01
export LTV_MAX=0.0
export GAMMA=5.0
export SUMMARY_JSON_PATH="$OUT_DIR/p6_option1_e1.json"

echo "=== Option-1 E1_2L run starting $(date) ===" | tee "$OUT_DIR/p6_option1_e1_stdout.log"
julia "$REPO_ROOT/src/vfi_solver_v4.jl" 2>&1 | tee -a "$OUT_DIR/p6_option1_e1_stdout.log"
echo "=== Done $(date) ===" | tee -a "$OUT_DIR/p6_option1_e1_stdout.log"
