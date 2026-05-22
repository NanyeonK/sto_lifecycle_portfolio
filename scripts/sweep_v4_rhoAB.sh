#!/usr/bin/env bash
# P1 sensitivity: rho_AB sweep for v4 Option 1.
# Tests whether CEV and hedge activation (mean_xB > 0) collapse as rho_AB → 1.
# At rho_AB → 1: cross-location returns perfectly correlated → no diversification
# benefit → pre-holding x_B should lose value → h1 hedge mechanism weakens.
#
# Run on server1 AFTER Step 6 (baseline runs) passes H2:
#   bash scripts/sweep_v4_rhoAB.sh
#
# Outputs: output/diagnostics/p6_v4_rhoAB/<REGIME>_rhoAB<tag>.json
#          output/diagnostics/p6_v4_rhoAB/summary.md

set -euo pipefail
JULIA=${JULIA:-julia}
SOLVER=src/vfi_solver_v4.jl
OUTDIR=output/diagnostics/p6_v4_rhoAB
mkdir -p "$OUTDIR"

RHO_AB_VALS="0.00 0.25 0.50 0.75 0.95"
COMMON_ENV="N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 ASSET_GRID_SIZE=7 GH_NODES=3 \
            TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
            P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
            GAMMA=5.0 BETA=0.96 LTV_MAX=0.0"

for rho_ab in $RHO_AB_VALS; do
    tag=$(echo "$rho_ab" | tr '.' 'p')
    for regime in E1_2L E2_2L; do
        outfile="$OUTDIR/${regime}_rhoAB${tag}.json"
        if [ -f "$outfile" ]; then
            echo "=== SKIP (exists): rho_AB=$rho_ab regime=$regime ==="
            continue
        fi
        echo "=== rho_AB=$rho_ab regime=$regime: $(date) ==="
        env $COMMON_ENV \
            REGIME="$regime" \
            RHO_AB="$rho_ab" \
            SUMMARY_JSON_PATH="$outfile" \
            "$JULIA" "$SOLVER" 2>&1 | tee "${outfile%.json}_stdout.log"
        echo "    done → $outfile"
    done
done

echo ""
echo "=== Computing CEV rhoAB sweep table ==="
"$JULIA" scripts/compute_cev_v4.jl rhoAB "$OUTDIR" \
    | tee "$OUTDIR/summary.md"
echo ""
echo "Done. Results in $OUTDIR/summary.md"
