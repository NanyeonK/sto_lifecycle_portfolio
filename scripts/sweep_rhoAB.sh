#!/usr/bin/env bash
# P1 sensitivity: rho_AB sweep using v4 solver (Option 1 6D state).
# Sweeps rho_AB ∈ {0, 0.25, 0.50, 0.75, 0.95} for E1_2L and E2_2L.
# At rho_AB → 1, the cross-location hedge should collapse.
#
# Prerequisites: v4 baselines (run_option1_e1.sh + run_option1_e2.sh) should
# have confirmed H1+H2 before this sweep is meaningful.
#
# Run on server1:
#   bash scripts/sweep_rhoAB.sh
#
# Outputs: output/diagnostics/p7_rhoAB_v4/<REGIME>_rhoAB<tag>.json
#          output/diagnostics/p7_rhoAB_v4/summary.md

set -euo pipefail
JULIA=${JULIA:-julia}
SOLVER=src/vfi_solver_v4.jl
OUTDIR=output/diagnostics/p7_rhoAB_v4
mkdir -p "$OUTDIR"

RHO_AB_VALS="0.00 0.25 0.50 0.75 0.95"

for rho_ab in $RHO_AB_VALS; do
    tag=$(echo "$rho_ab" | tr '.' 'p')
    for regime in E1_2L E2_2L; do
        outfile="$OUTDIR/${regime}_rhoAB${tag}.json"
        [ -f "$outfile" ] && { echo "SKIP: $outfile"; continue; }
        echo "=== rho_AB=$rho_ab regime=$regime ==="
        REGIME="$regime" \
        RHO_AB="$rho_ab" \
        N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=2.0 \
        ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
        GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
        SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 RHO=0.05 M=0.01 \
        TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
        P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
        LTV_MAX=0.0 \
        SUMMARY_JSON_PATH="$outfile" \
            "$JULIA" "$SOLVER" 2>&1 | tee "${outfile%.json}_stdout.log"
    done
done

echo ""
echo "=== CEV sweep table ==="
"$JULIA" scripts/compute_cev_sweep.jl "$OUTDIR" rhoAB | tee "$OUTDIR/summary.md"
echo "Done. $OUTDIR/summary.md"
