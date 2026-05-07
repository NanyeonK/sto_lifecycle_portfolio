#!/usr/bin/env bash
# run_option1_e2.sh — E2_2L Option 1 run (6D state, v4 solver, proper tau_buy)
# Run on server1 from repo root: bash scripts/run_option1_e2.sh
#
# Expected wall time: ~2-3h single thread (N_W=15, N_Z=5, N_X_PREV=3, GH=3)
# Output: output/diagnostics/p6_option1_e2.json  + stdout log
#
# Hypothesis 1: mean_xB_t1_xprev0_ellA > 0 (cross-location hedge activates)
# Hypothesis 2: CEV(E2_2L_v4 vs E1_2L_v4) > 4.255% (Option 3 baseline)
# Compute H3 (hedge channel CEV) by comparing V values between E2 and E1 outputs.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="$REPO_ROOT/output/diagnostics"
mkdir -p "$OUT_DIR"

LOG="$OUT_DIR/p6_option1_e2_stdout.log"
JSON="$OUT_DIR/p6_option1_e2.json"

echo "=== Option 1 E2_2L baseline (v4 solver, proper tau_buy) ===" | tee "$LOG"
date | tee -a "$LOG"

REGIME=E2_2L \
N_W=15 \
N_Z=5 \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5.0 \
BETA=0.96 \
SUMMARY_JSON_PATH="$JSON" \
  julia "$REPO_ROOT/src/vfi_solver_v4.jl" 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== E2_2L run complete ===" | tee -a "$LOG"
date | tee -a "$LOG"
echo "Summary JSON: $JSON" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "Next: check mean_xB_t1_xprev0_ellA in $JSON for hedge activation." | tee -a "$LOG"
echo "Then compute CEV(E2_2L_v4 vs E1_2L_v4) from V_t1_midpoint values." | tee -a "$LOG"
