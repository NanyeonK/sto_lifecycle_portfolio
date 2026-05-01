#!/usr/bin/env bash
# p4_full_txcost.sh — full round-trip transaction cost baseline (Round 4 P0 item 3).
#
# Runs E1_2L and E2_2L with full round-trip costs:
#   tau_sell = 0.06  (NAR sell-side ~6%)
#   tau_buy  = 0.025 (closing costs ~2.5%; now active Phase 1 approximation)
#   Total E1_2L round-trip (owning → relocate → own): 8.5% of housing value
#
# Also runs the channel decomposition counterfactual (tau_sell=0 for E1_2L_NOTX).
#
# Run from repo root on server1:
#   bash scripts/p4_full_txcost.sh 2>&1 | tee output/logs/p4_full_txcost.log
#
# Grid: N_W=21, N_Z=7, ASSET_GRID_SIZE=9, X_GRID_SIZE=7, GH_NODES=3

set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p output/diagnostics output/logs

GRID_ENV="N_W=21 N_Z=7 ASSET_GRID_SIZE=9 X_GRID_SIZE=7 GH_NODES=3"

echo "=== p4_full_txcost.sh  $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "Full round-trip: tau_sell=0.06, tau_buy=0.025 (active)"
echo ""

# 1. E1_2L baseline with full round-trip costs
echo "--- E1_2L (tau_sell=0.06, tau_buy=0.025) ---"
env $GRID_ENV \
    REGIME=E1_2L TAU_SELL=0.06 TAU_BUY=0.025 \
    SUMMARY_JSON_PATH=output/diagnostics/p4_txcost_E1_2L.json \
    julia src/vfi_solver_v3.jl

# 2. E1_2L_NOTX counterfactual (tau_sell=0, tau_buy=0; for channel decomp)
echo "--- E1_2L_NOTX (tau_sell=0, tau_buy=0) ---"
env $GRID_ENV \
    REGIME=E1_2L TAU_SELL=0.0 TAU_BUY=0.0 \
    SUMMARY_JSON_PATH=output/diagnostics/p4_txcost_E1_2L_NOTX.json \
    julia src/vfi_solver_v3.jl

# 3. E2_2L baseline (tokens: no forced sale, no buying cost)
echo "--- E2_2L (no relocation cost; tokens portable) ---"
env $GRID_ENV \
    REGIME=E2_2L TAU_SELL=0.06 TAU_BUY=0.025 \
    SUMMARY_JSON_PATH=output/diagnostics/p4_txcost_E2_2L.json \
    julia src/vfi_solver_v3.jl

echo ""
echo "=== CEV decomposition (full round-trip) ==="
julia scripts/compute_cev.jl \
    output/diagnostics/p4_txcost_E1_2L.json \
    output/diagnostics/p4_txcost_E2_2L.json \
    "TOTAL: CEV(E2_2L vs E1_2L) with full tx cost (tau_sell=6%+tau_buy=2.5%)"

julia scripts/compute_cev.jl \
    output/diagnostics/p4_txcost_E1_2L.json \
    output/diagnostics/p4_txcost_E1_2L_NOTX.json \
    "avoided-tx channel: CEV(E1_2L_NOTX vs E1_2L)"

julia scripts/compute_cev.jl \
    output/diagnostics/p4_txcost_E1_2L_NOTX.json \
    output/diagnostics/p4_txcost_E2_2L.json \
    "maintained-hedge channel: CEV(E2_2L vs E1_2L_NOTX)"

echo "=== p4_full_txcost.sh DONE ==="
echo "Summary JSONs: output/diagnostics/p4_txcost_*.json"
echo "Write results to: output/diagnostics/p4_full_txcost.md"
