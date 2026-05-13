#!/usr/bin/env bash
# P1 robustness: mortgage activation sweep using v4 solver (Option 1 6D state).
# Sweeps LTV_MAX ∈ {0.0, 0.50, 0.80} for E1_2L and E2_2L.
#
# Economic question: does the hedge channel (pre-buying x_B at ell=A)
# survive when households can leverage via mortgage?
# Expected pattern from v2: mortgage shrinks CEV(E2 vs E1) by ~37% because
# E1 households lever up via mortgage instead of needing continuous-x.
# v4 prediction: mortgage reduces continuous-x channel but the pre-buying
# hedge channel (uniquely unavailable in E1) should be more robust.
#
# Run on server1:
#   bash scripts/sweep_mortgage.sh
#
# Outputs: output/diagnostics/p7_mortgage_v4/<REGIME>_ltv<tag>.json
#          output/diagnostics/p7_mortgage_v4/summary.md

set -euo pipefail
JULIA=${JULIA:-julia}
SOLVER=src/vfi_solver_v4.jl
OUTDIR=output/diagnostics/p7_mortgage_v4
mkdir -p "$OUTDIR"

LTV_VALS="0.00 0.50 0.80"

for ltv in $LTV_VALS; do
    tag=$(echo "$ltv" | tr '.' 'p')
    for regime in E1_2L E2_2L; do
        outfile="$OUTDIR/${regime}_ltv${tag}.json"
        [ -f "$outfile" ] && { echo "SKIP: $outfile"; continue; }
        echo "=== LTV_MAX=$ltv regime=$regime ==="
        REGIME="$regime" \
        LTV_MAX="$ltv" \
        N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=2.0 \
        ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
        GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
        SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 RHO=0.05 M=0.01 \
        TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
        P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
        RHO_AB=0.50 \
        R_MORT_PREMIUM=0.005 \
        SUMMARY_JSON_PATH="$outfile" \
            "$JULIA" "$SOLVER" 2>&1 | tee "${outfile%.json}_stdout.log"
    done
done

echo ""
echo "=== CEV mortgage table ==="
"$JULIA" scripts/compute_cev_sweep.jl "$OUTDIR" mortgage | tee "$OUTDIR/summary.md"
echo "Done. $OUTDIR/summary.md"
