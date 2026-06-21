#!/usr/bin/env bash
# Run E2_2L baseline with v4 solver (Option 1: 6D state, proper tau_buy).
# Uses same grid settings as run_option1_e1.sh for direct CEV comparison.
# Expected wall time: ~60-90 min on server1 (single thread).
#
# Key test after run:
#   mean_xB_t1_ellA > 0  → hedge mechanism activated (H1 hypothesis)
#   CEV(E2_2L_v4 vs E1_2L_v4) > 4.255%  → H2 hypothesis
#   CEV(E2_2L_v4 vs E2_2L_v3) ≈ 0.5-1.5%  → H3 hypothesis (hedge channel)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$SCRIPT_DIR")"
mkdir -p "$REPO/output/diagnostics"

SUMMARY_JSON="$REPO/output/diagnostics/p6_option1_e2.json"
LOG_FILE="$REPO/output/diagnostics/p6_option1_e2_stdout.log"

echo "=== v4 E2_2L Option 1 run: $(date) ===" | tee "$LOG_FILE"

REGIME=E2_2L \
N_W=15 N_Z=5 \
N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 \
GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
RHO=0.05 M=0.01 SIGMA_H=0.115 SIGMA_DIV=0.10 \
SUMMARY_JSON_PATH="$SUMMARY_JSON" \
julia "$REPO/src/vfi_solver_v4.jl" 2>&1 | tee -a "$LOG_FILE"

echo "=== Done: $(date) ===" | tee -a "$LOG_FILE"
echo "Summary: $SUMMARY_JSON"
