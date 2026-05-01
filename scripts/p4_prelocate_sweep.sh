#!/usr/bin/env bash
# p4_prelocate_sweep.sh — p_relocate sensitivity sweep for Round 4 referee requirement.
#
# Sweeps P_RELOCATE_WORKING ∈ {0, 0.02, 0.06, 0.12} for E1_2L and E2_2L.
# Referee prediction: at p_relocate=0, cross-location holding must collapse
# (no relocation risk → no mobility-hedge motive).
#
# Run from repo root on server1:
#   bash scripts/p4_prelocate_sweep.sh 2>&1 | tee output/logs/p4_prelocate_sweep.log
#
# Grid: N_W=21, N_Z=7, ASSET_GRID_SIZE=9, X_GRID_SIZE=7, GH_NODES=3

set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p output/diagnostics output/logs

P_RELOC_VALUES=(0.0 0.02 0.06 0.12)
GRID_ENV="N_W=21 N_Z=7 ASSET_GRID_SIZE=9 X_GRID_SIZE=7 GH_NODES=3"

echo "=== p4_prelocate_sweep.sh  $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "Sweep: P_RELOCATE_WORKING in {${P_RELOC_VALUES[*]}}"
echo ""

for p_reloc in "${P_RELOC_VALUES[@]}"; do
    tag="preloc$(echo $p_reloc | tr '.' 'p')"   # e.g., preloc0p06
    echo "--- P_RELOCATE_WORKING = $p_reloc ---"

    for regime in E1_2L E2_2L; do
        json="output/diagnostics/p4_${tag}_${regime}.json"
        echo "  REGIME=$regime  P_RELOCATE_WORKING=$p_reloc  -> $json"
        env $GRID_ENV \
            REGIME="$regime" \
            P_RELOCATE_WORKING="$p_reloc" \
            SUMMARY_JSON_PATH="$json" \
            julia src/vfi_solver_v3.jl
    done

    # CEV at this p_relocate
    julia scripts/compute_cev.jl \
        "output/diagnostics/p4_${tag}_E1_2L.json" \
        "output/diagnostics/p4_${tag}_E2_2L.json" \
        "p_relocate_working=${p_reloc}"
done

echo "=== p4_prelocate_sweep.sh DONE ==="
echo "Summary JSONs: output/diagnostics/p4_preloc_*.json"
echo "Write results to: output/diagnostics/p4_prelocate_sweep.md"
