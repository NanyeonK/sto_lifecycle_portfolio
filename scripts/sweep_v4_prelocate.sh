#!/usr/bin/env bash
# P1 sensitivity: p_relocate sweep for v4 Option 1.
# Tests whether hedge activation scales with relocation frequency.
# At p_relocate=0: no relocation risk → pre-buying x_B has no value → mean_xB→0.
# At p_relocate=0.30: high mobility → strong pre-buy incentive → mean_xB↑.
#
# Run on server1 AFTER Step 6 (baseline runs) passes H2:
#   bash scripts/sweep_v4_prelocate.sh
#
# Outputs: output/diagnostics/p6_v4_prelocate/<REGIME>_preloc<tag>.json
#          output/diagnostics/p6_v4_prelocate/summary.md

set -euo pipefail
JULIA=${JULIA:-julia}
SOLVER=src/vfi_solver_v4.jl
OUTDIR=output/diagnostics/p6_v4_prelocate
mkdir -p "$OUTDIR"

PRELOC_VALS="0.00 0.06 0.12 0.30"
COMMON_ENV="N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 ASSET_GRID_SIZE=7 GH_NODES=3 \
            TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
            RHO_AB=0.50 GAMMA=5.0 BETA=0.96 LTV_MAX=0.0"

for preloc in $PRELOC_VALS; do
    tag=$(echo "$preloc" | tr '.' 'p')
    for regime in E1_2L E2_2L; do
        outfile="$OUTDIR/${regime}_preloc${tag}.json"
        if [ -f "$outfile" ]; then
            echo "=== SKIP (exists): p_reloc=$preloc regime=$regime ==="
            continue
        fi
        echo "=== p_relocate=$preloc regime=$regime: $(date) ==="
        env $COMMON_ENV \
            REGIME="$regime" \
            P_RELOCATE_WORKING="$preloc" \
            P_RELOCATE_RETIRED="$(awk "BEGIN{printf \"%.2f\", $preloc/3}")" \
            SUMMARY_JSON_PATH="$outfile" \
            "$JULIA" "$SOLVER" 2>&1 | tee "${outfile%.json}_stdout.log"
        echo "    done → $outfile"
    done
done

echo ""
echo "=== Computing CEV p_relocate sweep table ==="
"$JULIA" scripts/compute_cev_v4.jl prelocate "$OUTDIR" \
    | tee "$OUTDIR/summary.md"
echo ""
echo "Done. Results in $OUTDIR/summary.md"
