#!/usr/bin/env bash
# run_option1_e1_notx.sh — E1_2L counterfactual: NO transaction costs.
#
# Purpose: channel decomposition counterfactual for plot_channel_decomp.py.
#   CEV(E1_NOTX vs E1) = avoided-transaction-cost channel (tau_sell + tau_buy burden in E1).
#
# All params identical to run_option1_e1.sh EXCEPT:
#   TAU_SELL=0.0   (no selling cost on relocation)
#   TAU_BUY=0.0    (no buying cost on purchase)
#   TAU_TOKEN=0.0  (no transfer cost)
#
# Output: output/diagnostics/p6_option1_e1_notx.json
#
# Usage: bash scripts/run_option1_e1_notx.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$REPO_ROOT/output/diagnostics"

LOG="$REPO_ROOT/output/diagnostics/p6_option1_e1_notx_stdout.log"
JSON="$REPO_ROOT/output/diagnostics/p6_option1_e1_notx.json"

echo "=== v4 E1_2L_NOTX counterfactual run ==="
echo "  tau_sell=0, tau_buy=0, tau_token=0"
echo "  log  -> $LOG"
echo "  json -> $JSON"

REGIME=E1_2L \
N_W=15 N_Z=5 \
N_X_PREV=3 X_PREV_MAX=2.0 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 \
RHO=0.05 M=0.01 \
TAU_SELL=0.0 TAU_BUY=0.0 TAU_TOKEN=0.0 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$JSON" \
julia "$REPO_ROOT/src/vfi_solver_v4.jl" 2>&1 | tee "$LOG"

echo "=== E1_2L_NOTX done. Summary at $JSON ==="
