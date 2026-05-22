#!/usr/bin/env bash
# P1 sensitivity: LTV (mortgage) sweep for v4 Option 1.
# Mortgage availability substitutes for fractional ownership as a leveraging tool.
# Tests whether CEV shrinks as LTV increases (YZ/Cocco robustness check).
#
# Run on server1 AFTER Step 6 (baseline runs):
#   bash scripts/sweep_v4_mortgage.sh
#
# Outputs: output/diagnostics/p6_v4_mortgage/<REGIME>_ltv<tag>.json
#          output/diagnostics/p6_v4_mortgage/summary.md

set -euo pipefail
JULIA=${JULIA:-julia}
SOLVER=src/vfi_solver_v4.jl
OUTDIR=output/diagnostics/p6_v4_mortgage
mkdir -p "$OUTDIR"

LTV_VALS="0.0 0.5 0.8"
COMMON_ENV="N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 ASSET_GRID_SIZE=7 GH_NODES=3 \
            TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
            P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
            RHO_AB=0.50 GAMMA=5.0 BETA=0.96 R_MORT_PREMIUM=0.005"

for ltv in $LTV_VALS; do
    tag=$(echo "$ltv" | tr '.' 'p')
    for regime in E1_2L E2_2L; do
        outfile="$OUTDIR/${regime}_ltv${tag}.json"
        if [ -f "$outfile" ]; then
            echo "=== SKIP (exists): ltv_max=$ltv regime=$regime ==="
            continue
        fi
        echo "=== ltv_max=$ltv regime=$regime: $(date) ==="
        env $COMMON_ENV \
            REGIME="$regime" \
            LTV_MAX="$ltv" \
            SUMMARY_JSON_PATH="$outfile" \
            "$JULIA" "$SOLVER" 2>&1 | tee "${outfile%.json}_stdout.log"
        echo "    done → $outfile"
    done
done

echo ""
echo "=== Computing CEV mortgage sweep table ==="
"$JULIA" scripts/compute_cev_v4.jl mortgage "$OUTDIR" \
    | tee "$OUTDIR/summary.md"
echo ""
echo "Done. Results in $OUTDIR/summary.md"
