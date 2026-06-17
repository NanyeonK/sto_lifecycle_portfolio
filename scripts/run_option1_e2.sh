#!/usr/bin/env bash
# Run E2_2L baseline under v4 (6D state, proper tau_buy).
# Spec: handoff/tau_buy_option1_spec.md
#
# Default grids: N_W=15, N_Z=5, N_X_PREV=3, ASSET_GRID_SIZE=9, GH_NODES=3
# Expected wall time: ~50 min single thread on server1
#
# Usage:
#   bash scripts/run_option1_e2.sh
#   N_X_PREV=5 bash scripts/run_option1_e2.sh

set -euo pipefail

OUT_DIR="output/diagnostics"
mkdir -p "$OUT_DIR"

REGIME=E2_2L \
N_W="${N_W:-15}" \
N_Z="${N_Z:-5}" \
N_X_PREV="${N_X_PREV:-3}" \
X_PREV_MAX="1.0" \
ASSET_GRID_SIZE="${ASSET_GRID_SIZE:-9}" \
GH_NODES="3" \
GAMMA="5.0" \
RHO="0.05" \
M="0.01" \
TAU_SELL="0.06" \
TAU_BUY="0.025" \
TAU_TOKEN="0.01" \
P_RELOCATE_WORKING="0.06" \
P_RELOCATE_RETIRED="0.02" \
RHO_AB="0.50" \
SIGMA_DIV="0.10" \
LTV_MAX="0.0" \
SUMMARY_JSON_PATH="$OUT_DIR/p6_option1_e2.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUT_DIR/p6_option1_e2_stdout.log"

echo "Done. Summary: $OUT_DIR/p6_option1_e2.json"
