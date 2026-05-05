#!/usr/bin/env bash
# P1 sensitivity: transaction-cost sweep using v4 solver.
# Varies (tau_sell, tau_buy) to isolate how much of the CEV comes from
# avoided round-trip costs vs the continuous-x / hedge channel.
#
# v4 note: tau_buy is applied natively via per-period delta (no APPLY_TAU_BUY flag).
# Setting TAU_BUY=0 completely removes the buying-cost channel from both regimes.
#
# Scenarios:
#   notx   — tau_sell=0,    tau_buy=0       (no tx costs; pure continuous-x)
#   sell6  — tau_sell=6%,   tau_buy=0       (selling friction only)
#   rt8p5  — tau_sell=6%,   tau_buy=2.5%   (8.5% round-trip; NAR baseline)
#   rt10   — tau_sell=6%,   tau_buy=4.0%   (10% round-trip)
#   rt12   — tau_sell=6%,   tau_buy=6.0%   (12% round-trip; high-cost market)
#
# Run on server1:
#   bash scripts/sweep_txcost.sh
#
# Outputs: output/diagnostics/p7_txcost_v4/<REGIME>_<tag>.json
#          output/diagnostics/p7_txcost_v4/summary.md

set -euo pipefail
JULIA=${JULIA:-julia}
SOLVER=src/vfi_solver_v4.jl
OUTDIR=output/diagnostics/p7_txcost_v4
mkdir -p "$OUTDIR"

COMMON="N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=2.0
        ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3
        GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04
        SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 RHO=0.05 M=0.01
        TAU_TOKEN=0.01 P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02
        RHO_AB=0.50 LTV_MAX=0.0"

run_regime () {
    local regime="$1"
    local tag="$2"
    local ts="$3"
    local tb="$4"
    local outfile="$OUTDIR/${regime}_${tag}.json"
    [ -f "$outfile" ] && { echo "SKIP: $outfile"; return; }
    echo "=== regime=$regime tag=$tag tau_sell=$ts tau_buy=$tb ==="
    env $COMMON \
        REGIME="$regime" TAU_SELL="$ts" TAU_BUY="$tb" \
        SUMMARY_JSON_PATH="$outfile" \
        "$JULIA" "$SOLVER" 2>&1 | tee "${outfile%.json}_stdout.log"
}

for regime in E1_2L E2_2L; do
    run_regime "$regime" "notx"   0.00  0.000
    run_regime "$regime" "sell6"  0.06  0.000
    run_regime "$regime" "rt8p5"  0.06  0.025
    run_regime "$regime" "rt10"   0.06  0.040
    run_regime "$regime" "rt12"   0.06  0.060
done

echo ""
echo "=== CEV sweep table ==="
"$JULIA" scripts/compute_cev_sweep.jl "$OUTDIR" txcost | tee "$OUTDIR/summary.md"
echo "Done. $OUTDIR/summary.md"
