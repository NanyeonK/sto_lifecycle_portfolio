#!/usr/bin/env bash
# Run E1_2L baseline with v4 solver (Option 1: full x_prev state)
# Usage: bash scripts/run_option1_e1.sh
# Server1 expected wall time: ~45 min (single thread)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$REPO_ROOT/output/diagnostics"
mkdir -p "$OUT"

export REGIME=E1_2L
export N_W=15
export N_Z=5
export N_X_PREV=3
export X_PREV_MAX=1.5
export ASSET_GRID_SIZE=9
export X_GRID_SIZE=5
export GH_NODES=3
export TAU_SELL=0.06
export TAU_BUY=0.025
export TAU_TOKEN=0.01
export P_RELOCATE_WORKING=0.06
export P_RELOCATE_RETIRED=0.02
export RHO_AB=0.50
export LTV_MAX=0.0
export SUMMARY_JSON_PATH="$OUT/p6_option1_e1.json"

echo "=== v4 E1_2L baseline (Option 1) ==="
julia "$REPO_ROOT/src/vfi_solver_v4.jl" 2>&1 | tee "$OUT/p6_option1_e1_stdout.log"
echo "Done. JSON: $SUMMARY_JSON_PATH"
