#!/usr/bin/env bash
# run_p4_full_txcost.sh — Round 4 P0: full round-trip transaction cost comparison
#
# Approximation for tau_buy (Phase 2 state-extension deferred):
#   E1_2L uses TAU_SELL = tau_sell + tau_buy = 0.06 + 0.025 = 0.085 (round-trip 8.5%)
#   E2_2L tokens are portable: TAU_SELL = 0.0 (no forced sale on relocation)
#   E1_2L_NOTX (counterfactual zero tx): TAU_SELL = 0.0 (for channel decomp comparison)
#
# Done artifact: output/diagnostics/p4_full_txcost.md
#
# Usage (server1, from repo root):
#   bash scripts/run_p4_full_txcost.sh
#
# Requires: julia, JSON3 package, src/vfi_solver_v3.jl
# Estimated wall time: ~30 min per regime at default grids, single thread.
# Run in separate tmux panes or with GNU parallel for parallel execution.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIAG_DIR="${REPO_ROOT}/output/diagnostics"
LOG_DIR="${REPO_ROOT}/output/logs"
SOLVER="${REPO_ROOT}/src/vfi_solver_v3.jl"

mkdir -p "${DIAG_DIR}" "${LOG_DIR}"

# Default calibration: full grids, single-thread
COMMON_VARS=(
    "JULIA_NUM_THREADS=1"
    "BETA=0.96"
    "GAMMA=5.0"
    "RF=1.02"
    "EQUITY_PREMIUM=0.04"
    "RHO=0.05"
    "M=0.01"
    "SIGMA_H=0.115"
    "SIGMA_DIV=0.10"
    "RHO_AB=0.50"
    "P_RELOCATE_WORKING=0.06"
    "P_RELOCATE_RETIRED=0.02"
    "TAU_TOKEN=0.01"
    "LTV_MAX=0.0"
)

run_regime() {
    local label="$1"
    local regime="$2"
    local tau_sell="$3"
    local out_json="${DIAG_DIR}/p4_full_txcost_${label}.json"
    local out_log="${LOG_DIR}/p4_full_txcost_${label}_stdout.log"

    echo "[$(date -u +%H:%M:%S)] Starting ${label} (REGIME=${regime}, TAU_SELL=${tau_sell})"
    env "${COMMON_VARS[@]}" \
        REGIME="${regime}" \
        TAU_SELL="${tau_sell}" \
        SUMMARY_JSON_PATH="${out_json}" \
        julia "${SOLVER}" 2>&1 | tee "${out_log}"
    echo "[$(date -u +%H:%M:%S)] Done ${label} → ${out_json}"
}

# ── Regime runs ───────────────────────────────────────────────────────────────
# 1. E1_2L baseline: sell 6% + buy 2.5% = 8.5% round-trip
run_regime "E1_2L_fulltx"   "E1_2L"  "0.085"

# 2. E2_2L: tokens portable, zero forced-sale cost
run_regime "E2_2L_notx"     "E2_2L"  "0.0"

# 3. E1_2L counterfactual zero tx (for channel decomposition comparison)
run_regime "E1_2L_notx"     "E1_2L"  "0.0"

echo ""
echo "=== All p4_full_txcost regimes complete ==="
echo "Compute CEV and write output/diagnostics/p4_full_txcost.md manually"
echo "using V values from the three JSON files above."
echo ""
echo "CEV formula (representative midpoint):"
echo "  V_E1_2L_fulltx = V from p4_full_txcost_E1_2L_fulltx.json"
echo "  V_E2_2L_notx   = V from p4_full_txcost_E2_2L_notx.json"
echo "  V_E1_2L_notx   = V from p4_full_txcost_E1_2L_notx.json"
echo ""
echo "  CEV_total          = (V_E2_2L_notx / V_E1_2L_fulltx)^(1/(1-gamma)) - 1"
echo "  CEV_avoided_tx     = (V_E1_2L_notx / V_E1_2L_fulltx)^(1/(1-gamma)) - 1"
echo "  CEV_hedge_channel  = (V_E2_2L_notx / V_E1_2L_notx)^(1/(1-gamma)) - 1"
echo "  cross_term         = CEV_total - CEV_avoided_tx - CEV_hedge_channel"
