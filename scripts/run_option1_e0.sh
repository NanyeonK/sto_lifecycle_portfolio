#!/usr/bin/env bash
# run_option1_e0.sh — E0 (rent-only) baseline with v4 solver.
#
# Purpose: E0 is the lowest benchmark regime (pure renter, no housing asset).
#   CEV(E1 vs E0) = value of binary homeownership (Cocco 2005 / YZ 2005 territory).
#   CEV(E2 vs E0) = total lifetime value of tokenized housing.
#   V_E0 is also needed for Figure 4 (V-slice three-regime comparison).
#
# Grid: same as e1/e2 runs for clean comparison.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$REPO_ROOT/output/diagnostics"

LOG="$REPO_ROOT/output/diagnostics/p6_option1_e0_stdout.log"
JSON="$REPO_ROOT/output/diagnostics/p6_option1_e0.json"

echo "=== v4 E0 (rent-only) baseline run ==="
echo "  log  -> $LOG"
echo "  json -> $JSON"

REGIME=E0 \
N_W=15 N_Z=5 \
N_X_PREV=3 X_PREV_MAX=2.0 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 \
RHO=0.05 M=0.01 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$JSON" \
julia "$REPO_ROOT/src/vfi_solver_v4.jl" 2>&1 | tee "$LOG"

echo "=== E0 done. Summary at $JSON ==="
