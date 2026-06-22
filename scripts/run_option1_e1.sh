#!/usr/bin/env bash
# run_option1_e1.sh — E1_2L baseline at v4 settings (Option 1 full state extension)
# Run on server1:  bash scripts/run_option1_e1.sh
# Output: output/diagnostics/p6_option1_e1.json + p6_option1_e1_stdout.log
#
# Grid settings (offset the 9x state factor from N_X_PREV=3):
#   N_W=15 (down from 21), N_Z=5 (down from 7)
#   N_X_PREV=3  →  x_prev_grid = {0.0, 0.5, 1.0}
#   Net compute: ~9 * (15*5)/(21*7) ≈ 4.6x v3 per regime; ~2-3h wall expected.

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

LOG="$OUTDIR/p6_option1_e1_stdout.log"
JSON="$OUTDIR/p6_option1_e1.json"

echo "=== Option 1 E1_2L run starting at $(date) ===" | tee "$LOG"

REGIME=E1_2L \
N_W=15 N_Z=5 \
N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
APPLY_TAU_BUY=0 \
RHO=0.05 M=0.01 \
SIGMA_H=0.115 SIGMA_DIV=0.10 \
RHO_AB=0.50 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$JSON" \
julia src/vfi_solver_v4.jl 2>&1 | tee -a "$LOG"

echo "=== E1_2L run complete at $(date) ===" | tee -a "$LOG"
echo "JSON summary: $JSON"
