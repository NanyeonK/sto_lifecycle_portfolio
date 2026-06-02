#!/usr/bin/env bash
# run_option1_e2_notau.sh — E2_2L counterfactual: NO tau_buy or tau_token.
#
# Purpose: channel decomposition counterfactual for plot_channel_decomp.py.
#   CEV(E2_NOTAU vs E1_NOTX) = continuous-x rent-saving channel (Liu 2021 territory).
#   CEV(E2      vs E2_NOTAU) = pre-buying hedge channel (uniquely tokens-enabled).
#
# All params identical to run_option1_e2.sh EXCEPT:
#   TAU_BUY=0.0    (no per-period buying cost on token increments)
#   TAU_TOKEN=0.0  (no per-period selling cost on token reductions)
#   TAU_SELL=0.06  (tau_sell retained: forced relocation sale still costs 6%,
#                   isolating the *pre-buying hedge* benefit from tau_buy portability)
#
# Output: output/diagnostics/p6_option1_e2_notau.json
#
# Usage: bash scripts/run_option1_e2_notau.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$REPO_ROOT/output/diagnostics"

LOG="$REPO_ROOT/output/diagnostics/p6_option1_e2_notau_stdout.log"
JSON="$REPO_ROOT/output/diagnostics/p6_option1_e2_notau.json"

echo "=== v4 E2_2L_NOTAU counterfactual run ==="
echo "  tau_buy=0, tau_token=0 (tau_sell retained at 0.06 for E1 forced-sale baseline parity)"
echo "  log  -> $LOG"
echo "  json -> $JSON"

REGIME=E2_2L \
N_W=15 N_Z=5 \
N_X_PREV=3 X_PREV_MAX=2.0 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 \
RHO=0.05 M=0.01 \
TAU_SELL=0.06 TAU_BUY=0.0 TAU_TOKEN=0.0 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$JSON" \
julia "$REPO_ROOT/src/vfi_solver_v4.jl" 2>&1 | tee "$LOG"

echo "=== E2_2L_NOTAU done. Summary at $JSON ==="
