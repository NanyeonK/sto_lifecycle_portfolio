#!/usr/bin/env bash
# run_p4_rhoAB_sweep.sh — Round 4 P1: cross-location correlation sensitivity
#
# Hypothesis: hedge channel collapses at rho_AB → 1 (perfect correlation ⟹
# x_B has same return as x_A; no diversification benefit from cross-location
# holding; household reverts to single-location allocation).
#
# Sweep: RHO_AB ∈ {0.00, 0.25, 0.50, 0.75, 0.95}
#   - Baseline is RHO_AB=0.50 (Case-Shiller MSA-pair midpoint).
#   - Both E1_2L and E2_2L are run at each value so that
#     CEV(E2_2L vs E1_2L) can be computed per rho_AB.
#
# Done artifact: output/diagnostics/p4_rhoAB_sweep.md
#
# Usage (server1, from repo root):
#   bash scripts/run_p4_rhoAB_sweep.sh
#   # or, for parallelism across RHO_AB values:
#   bash scripts/run_p4_rhoAB_sweep.sh &
#
# Estimated wall time: ~30 min per regime × 5 values × 2 regimes = 5 hours
# sequential. With 10 parallel workers, ~30 min total.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIAG_DIR="${REPO_ROOT}/output/diagnostics"
LOG_DIR="${REPO_ROOT}/output/logs"
SOLVER="${REPO_ROOT}/src/vfi_solver_v3.jl"

mkdir -p "${DIAG_DIR}" "${LOG_DIR}"

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
    "P_RELOCATE_WORKING=0.06"
    "P_RELOCATE_RETIRED=0.02"
    "TAU_SELL=0.06"
    "TAU_TOKEN=0.01"
    "LTV_MAX=0.0"
)

RHO_AB_VALUES=(0.00 0.25 0.50 0.75 0.95)
REGIMES=("E1_2L" "E2_2L")

for rho_ab in "${RHO_AB_VALUES[@]}"; do
    rho_tag="${rho_ab//./_}"   # e.g. 0.25 → 0_25
    for regime in "${REGIMES[@]}"; do
        label="p4_rhoAB_${rho_tag}_${regime}"
        out_json="${DIAG_DIR}/${label}.json"
        out_log="${LOG_DIR}/${label}_stdout.log"

        echo "[$(date -u +%H:%M:%S)] Starting ${label}"
        env "${COMMON_VARS[@]}" \
            REGIME="${regime}" \
            RHO_AB="${rho_ab}" \
            SUMMARY_JSON_PATH="${out_json}" \
            julia "${SOLVER}" 2>&1 | tee "${out_log}"
        echo "[$(date -u +%H:%M:%S)] Done ${label} → ${out_json}"
    done
done

echo ""
echo "=== rho_AB sweep complete ==="
echo "10 JSON files written to ${DIAG_DIR}/ (5 rho_AB values × 2 regimes)."
echo "Fill in output/diagnostics/p4_rhoAB_sweep.md with CEV values per rho_AB."
echo ""
echo "Key prediction: CEV(E2_2L vs E1_2L) shrinks monotonically as rho_AB → 1."
echo "At rho_AB = 0.95: cross-location tokens provide near-zero diversification;"
echo "maintained-hedge channel collapses toward avoided-tx channel only."
