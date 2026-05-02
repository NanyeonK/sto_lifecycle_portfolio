#!/usr/bin/env bash
# run_option1_e1.sh — E1_2L baseline at v4 Option 1 settings
# Usage: bash scripts/run_option1_e1.sh
#
# Coarse grid (N_X_PREV=3) recommended for first confirmation run.
# Expect ~2-3 hours wall on server1 single thread.

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

REGIME=E1_2L
SUMMARY_JSON="$OUTDIR/p6_option1_e1.json"
STDOUT_LOG="$OUTDIR/p6_option1_e1_stdout.log"

echo "=== v4 Option 1 E1_2L baseline ==="
echo "  REGIME         = $REGIME"
echo "  N_X_PREV       = ${N_X_PREV:-3}"
echo "  N_W            = ${N_W:-15}"
echo "  N_Z            = ${N_Z:-5}"
echo "  TAU_SELL       = ${TAU_SELL:-0.06}"
echo "  TAU_BUY        = ${TAU_BUY:-0.025}"
echo "  TAU_TOKEN      = ${TAU_TOKEN:-0.005}"
echo "  Output JSON    = $SUMMARY_JSON"
echo ""

REGIME=$REGIME \
  N_W=${N_W:-15} \
  N_Z=${N_Z:-5} \
  N_X_PREV=${N_X_PREV:-3} \
  X_PREV_MAX=${X_PREV_MAX:-1.0} \
  TAU_SELL=${TAU_SELL:-0.06} \
  TAU_BUY=${TAU_BUY:-0.025} \
  TAU_TOKEN=${TAU_TOKEN:-0.005} \
  APPLY_TAU_BUY=0 \
  SUMMARY_JSON_PATH="$SUMMARY_JSON" \
  julia src/vfi_solver_v4.jl 2>&1 | tee "$STDOUT_LOG"

echo ""
echo "Done. Summary at: $SUMMARY_JSON"
