#!/usr/bin/env bash
# Run E1_2L baseline for Path B Option 1 (v4 solver, 6D state).
# E1_2L = traditional binary homeownership; on relocation, forced sell (tau_sell)
# and forced buy (tau_buy) at new location.
# Expected runtime: ~30-60 min on server1 (single thread, coarse grid).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="$REPO_ROOT/output/diagnostics"
mkdir -p "$OUT_DIR"

echo "Starting E1_2L v4 baseline at $(date)"
SUMMARY_JSON_PATH="$OUT_DIR/p6_option1_e1.json" \
  REGIME=E1_2L \
  N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
  ASSET_GRID_SIZE=9 GH_NODES=3 \
  TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
  P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
  RHO_AB=0.5 LTV_MAX=0.0 \
  julia "$REPO_ROOT/src/vfi_solver_v4.jl" \
    2>&1 | tee "$OUT_DIR/p6_option1_e1_stdout.log"

echo "E1_2L done at $(date). JSON: $OUT_DIR/p6_option1_e1.json"
