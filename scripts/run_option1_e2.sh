#!/usr/bin/env bash
# Run E2_2L baseline for v4 (6D state / Option 1 full tx_cost)
# Spec: handoff/tau_buy_option1_spec.md
# Run on server1: bash scripts/run_option1_e2.sh
# Output: output/diagnostics/p6_option1_e2.json

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="$REPO_DIR/output/diagnostics"
mkdir -p "$OUT_DIR"

SUMMARY_JSON_PATH="$OUT_DIR/p6_option1_e2.json" \
REGIME="E2_2L" \
N_W=15 \
N_Z=5 \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 \
GH_NODES=3 \
GAMMA=5.0 \
BETA=0.96 \
RF=1.02 \
EQUITY_PREMIUM=0.04 \
SIGMA_H=0.115 \
SIGMA_DIV=0.10 \
G_H=0.016 \
RHO=0.05 \
M=0.01 \
RHO_AB=0.5 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.005 \
LTV_MAX=0.0 \
SAVE_PATH="$OUT_DIR/p6_option1_e2_result.jls" \
  julia "$REPO_DIR/src/vfi_solver_v4.jl" 2>&1 | tee "$OUT_DIR/p6_option1_e2_stdout.log"

echo "E2_2L v4 run complete. Summary: $OUT_DIR/p6_option1_e2.json"
