#!/usr/bin/env bash
# run_p4_prelocate_sweep.sh — Round 4 P1: relocation probability sensitivity
#
# Hypothesis: cross-location token holding collapses at p_relocate = 0.
# If the household never relocates, the token portability mechanism has no
# value (no relocation events to hedge against), and E2_2L reverts to E1_2L.
# CEV(E2_2L vs E1_2L) should shrink monotonically toward zero as p → 0.
#
# Sweep: P_RELOCATE_WORKING ∈ {0.00, 0.02, 0.06, 0.12}
#   - p=0.00: no relocation (should give CEV ≈ 0 from avoided-tx only)
#   - p=0.02: low mobility (1-in-50 annual)
#   - p=0.06: baseline (PSID mid-range working-age)
#   - p=0.12: high mobility (urban movers, PSID upper range)
#
#   P_RELOCATE_RETIRED fixed at 0.02 throughout (retirement mobility is low).
#
# Done artifact: output/diagnostics/p4_prelocate_sweep.md
#
# Usage (server1, from repo root):
#   bash scripts/run_p4_prelocate_sweep.sh
#
# Estimated wall time: ~30 min per regime × 4 values × 2 regimes = 4 hours
# sequential. With 8 parallel workers, ~30 min total.

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
    "RHO_AB=0.50"
    "P_RELOCATE_RETIRED=0.02"
    "TAU_SELL=0.06"
    "TAU_TOKEN=0.01"
    "LTV_MAX=0.0"
)

P_WORK_VALUES=(0.00 0.02 0.06 0.12)
REGIMES=("E1_2L" "E2_2L")

for p_work in "${P_WORK_VALUES[@]}"; do
    p_tag="${p_work//./_}"   # e.g. 0.06 → 0_06
    for regime in "${REGIMES[@]}"; do
        label="p4_prelocate_${p_tag}_${regime}"
        out_json="${DIAG_DIR}/${label}.json"
        out_log="${LOG_DIR}/${label}_stdout.log"

        echo "[$(date -u +%H:%M:%S)] Starting ${label}"
        env "${COMMON_VARS[@]}" \
            REGIME="${regime}" \
            P_RELOCATE_WORKING="${p_work}" \
            SUMMARY_JSON_PATH="${out_json}" \
            julia "${SOLVER}" 2>&1 | tee "${out_log}"
        echo "[$(date -u +%H:%M:%S)] Done ${label} → ${out_json}"
    done
done

echo ""
echo "=== p_relocate sweep complete ==="
echo "8 JSON files written to ${DIAG_DIR}/ (4 p_work values × 2 regimes)."
echo "Fill in output/diagnostics/p4_prelocate_sweep.md with CEV values."
echo ""
echo "Key predictions:"
echo "  - At p=0.00: CEV(E2_2L vs E1_2L) ≈ 0 (no relocation events to hedge)"
echo "  - At p=0.06 (baseline): CEV ≈ +4.23% (from channel decomposition run)"
echo "  - At p=0.12: CEV > baseline (more relocations amplify both channels)"
echo "  - mean_xB (cross-location token share) should collapse at p=0.00"
