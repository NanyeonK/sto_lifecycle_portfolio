#!/usr/bin/env bash
# run_option1_e2.sh — Run E2_2L baseline with v4 6D state (Option 1 tau_buy)
#
# Usage: bash scripts/run_option1_e2.sh
# Output: output/diagnostics/p6_option1_e2.json  +  p6_option1_e2_stdout.log
#
# Grid: N_W=15, N_Z=5, N_X_PREV=3 (matches e1 run for clean CEV comparison)
# Regime: E2_2L — continuous fractional tokens; tokens portable across relocation;
#         tau_buy charged on any positive delta in x_A or x_B each period.
#         Pre-holding x_B at ell=A reduces tau_buy on relocation to B (the hedge channel).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$REPO_ROOT/output/diagnostics"

LOG="$REPO_ROOT/output/diagnostics/p6_option1_e2_stdout.log"
JSON="$REPO_ROOT/output/diagnostics/p6_option1_e2.json"

echo "=== v4 E2_2L Option 1 run ==="
echo "  log  -> $LOG"
echo "  json -> $JSON"

REGIME=E2_2L \
N_W=15 N_Z=5 \
N_X_PREV=3 X_PREV_MAX=2.0 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 \
RHO=0.05 M=0.01 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$JSON" \
julia "$REPO_ROOT/src/vfi_solver_v4.jl" 2>&1 | tee "$LOG"

echo "=== E2_2L done. Summary at $JSON ==="
