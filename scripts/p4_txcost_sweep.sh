#!/usr/bin/env bash
# p4_txcost_sweep.sh — Round 4 MUST P0: full transaction-cost sweep
#
# Runs E1_2L and E2_2L under three tau configurations to quantify:
#   (a) avoided-tx channel (tau_sell + tau_buy impact on E1_2L)
#   (b) CEV(E2_2L vs E1_2L) under realistic round-trip costs
#
# Requires: Julia installed, vfi_solver_v3.jl in ../src/
# Run on server1:
#   bash scripts/p4_txcost_sweep.sh 2>&1 | tee output/logs/p4_txcost_sweep.log
#
# Results: JSON summaries written to output/diagnostics/
#          Combine manually with scripts/compute_cev_v3.jl
#
# Three scenarios:
#   (1) tx_baseline:  tau_sell=0.06, tau_buy=0.025  (NAR sell + typical closing costs)
#   (2) tx_roundtrip: tau_sell=0.085, tau_buy=0.025  (8.5% sell + 2.5% buy = ~11% round-trip)
#   (3) tx_none:      tau_sell=0.00, tau_buy=0.00    (no tx costs — counterfactual)
#
# CEV decomposition:
#   TOTAL tokenization gain:      CEV(E2_2L vs E1_2L | tx_baseline)
#   Avoided-tx channel:           CEV(E1_2L_tx_none vs E1_2L | tx_baseline)
#   Maintained-hedge channel:     CEV(E2_2L vs E1_2L_tx_none)
#   Round-trip sensitivity:       CEV(E2_2L vs E1_2L | tx_roundtrip)

set -euo pipefail

JULIA="${JULIA_BIN:-julia}"
SOLVER="$(dirname "$0")/../src/vfi_solver_v3.jl"
OUTDIR="$(dirname "$0")/../output/diagnostics"
mkdir -p "$OUTDIR"

# Grid settings — use default (full-grid) unless overridden
N_W="${N_W:-81}"
N_Z="${N_Z:-11}"
ASSET_GRID_SIZE="${ASSET_GRID_SIZE:-21}"
X_GRID_SIZE="${X_GRID_SIZE:-11}"
GH_NODES="${GH_NODES:-3}"

export N_W N_Z ASSET_GRID_SIZE X_GRID_SIZE GH_NODES

echo "=== p4_txcost_sweep.sh — $(date) ==="
echo "  N_W=$N_W  N_Z=$N_Z  ASSET=$ASSET_GRID_SIZE  X=$X_GRID_SIZE  GH=$GH_NODES"

# ── Scenario 1: tx_baseline (tau_sell=6%, tau_buy=2.5%) ──────────────────────

echo "--- [1/6] E1_2L tx_baseline ---"
REGIME=E1_2L TAU_SELL=0.06 TAU_BUY=0.025 \
  SUMMARY_JSON_PATH="$OUTDIR/p4_txcost_E1_2L_baseline.json" \
  "$JULIA" "$SOLVER" \
  2>&1 | tee "$OUTDIR/p4_txcost_E1_2L_baseline.log"

echo "--- [2/6] E2_2L tx_baseline ---"
REGIME=E2_2L TAU_SELL=0.06 TAU_BUY=0.025 \
  SUMMARY_JSON_PATH="$OUTDIR/p4_txcost_E2_2L_baseline.json" \
  "$JULIA" "$SOLVER" \
  2>&1 | tee "$OUTDIR/p4_txcost_E2_2L_baseline.log"

# ── Scenario 2: tx_roundtrip (tau_sell=8.5%, tau_buy=2.5%) ───────────────────
# Represents upper NAR range: 7-8% sell + 2-3% buy = ~10-11% round-trip.

echo "--- [3/6] E1_2L tx_roundtrip ---"
REGIME=E1_2L TAU_SELL=0.085 TAU_BUY=0.025 \
  SUMMARY_JSON_PATH="$OUTDIR/p4_txcost_E1_2L_roundtrip.json" \
  "$JULIA" "$SOLVER" \
  2>&1 | tee "$OUTDIR/p4_txcost_E1_2L_roundtrip.log"

echo "--- [4/6] E2_2L tx_roundtrip ---"
REGIME=E2_2L TAU_SELL=0.085 TAU_BUY=0.025 \
  SUMMARY_JSON_PATH="$OUTDIR/p4_txcost_E2_2L_roundtrip.json" \
  "$JULIA" "$SOLVER" \
  2>&1 | tee "$OUTDIR/p4_txcost_E2_2L_roundtrip.log"

# ── Scenario 3: tx_none — no-tx counterfactual (matches E1_2L_NOTX) ─────────

echo "--- [5/6] E1_2L tx_none ---"
REGIME=E1_2L TAU_SELL=0.00 TAU_BUY=0.00 \
  SUMMARY_JSON_PATH="$OUTDIR/p4_txcost_E1_2L_notx.json" \
  "$JULIA" "$SOLVER" \
  2>&1 | tee "$OUTDIR/p4_txcost_E1_2L_notx.log"

echo "--- [6/6] E2_2L tx_none ---"
REGIME=E2_2L TAU_SELL=0.00 TAU_BUY=0.00 \
  SUMMARY_JSON_PATH="$OUTDIR/p4_txcost_E2_2L_notx.json" \
  "$JULIA" "$SOLVER" \
  2>&1 | tee "$OUTDIR/p4_txcost_E2_2L_notx.log"

echo ""
echo "=== All runs complete. JSON summaries in $OUTDIR/ ==="
echo "    Next: run scripts/compute_cev_v3.jl to compile CEV decomposition table"
echo "    Expected output: output/diagnostics/p4_full_txcost.md"
