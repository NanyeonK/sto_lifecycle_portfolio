#!/usr/bin/env bash
# p4_rhoAB_sweep.sh — rho_AB sensitivity sweep for Round 4 referee requirement.
#
# Sweeps rho_AB ∈ {0, 0.25, 0.5, 0.75, 0.95} for E1_2L and E2_2L.
# Referee prediction: hedge channel collapses as rho_AB → 1
# (cross-location diversification disappears at perfect correlation).
#
# Run from repo root on server1:
#   bash scripts/p4_rhoAB_sweep.sh 2>&1 | tee output/logs/p4_rhoAB_sweep.log
#
# Grid: N_W=21, N_Z=7, ASSET_GRID_SIZE=9, X_GRID_SIZE=7, GH_NODES=3
# (matches previous full-grid baseline runs).

set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p output/diagnostics output/logs

RHO_AB_VALUES=(0.0 0.25 0.5 0.75 0.95)
GRID_ENV="N_W=21 N_Z=7 ASSET_GRID_SIZE=9 X_GRID_SIZE=7 GH_NODES=3"

echo "=== p4_rhoAB_sweep.sh  $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "Sweep: rho_AB in {${RHO_AB_VALUES[*]}}"
echo ""

for rho_ab in "${RHO_AB_VALUES[@]}"; do
    tag="rhoAB$(echo $rho_ab | tr '.' 'p')"   # e.g., rhoAB0p5
    echo "--- rho_AB = $rho_ab ---"

    for regime in E1_2L E2_2L; do
        json="output/diagnostics/p4_${tag}_${regime}.json"
        echo "  REGIME=$regime  rho_AB=$rho_ab  -> $json"
        env $GRID_ENV \
            REGIME="$regime" \
            RHO_AB="$rho_ab" \
            SUMMARY_JSON_PATH="$json" \
            julia src/vfi_solver_v3.jl
    done

    # CEV at this rho_AB
    julia scripts/compute_cev.jl \
        "output/diagnostics/p4_${tag}_E1_2L.json" \
        "output/diagnostics/p4_${tag}_E2_2L.json" \
        "rho_AB=${rho_ab}"
done

echo "=== p4_rhoAB_sweep.sh DONE ==="
echo "Summary JSONs: output/diagnostics/p4_rhoAB_*.json"
echo "Write results to: output/diagnostics/p4_rhoAB_sweep.md"
