#!/usr/bin/env bash
# Round 4 P0-3: tau_buy / round-trip transaction cost sweep.
# Tests E1_2L vs E2_2L under a range of round-trip cost scenarios:
#   - NOTX: tau_sell=0, tau_buy=0 (counterfactual, replicates channel-decomp E1_2L_NOTX)
#   - sell6:      tau_sell=6%, tau_buy=0    (sell-only, current baseline)
#   - rt8p5:      tau_sell=6%, tau_buy=2.5% (NAR sell + typical closing costs, 8.5% rt)
#   - rt10:       tau_sell=6%, tau_buy=4%   (upper NAR range, 10% rt)
#   - rt12:       tau_sell=6%, tau_buy=6%   (high-cost market, 12% rt)
#
# tau_buy is applied via the apply_tau_buy_at_reloc approximation (see solver comments):
# owner who relocates pays tau_sell on sell proceeds AND tau_buy deducted from reloc wealth,
# assuming they re-purchase 1 unit at the new location. E2_2L: no reloc costs (tokens portable).
#
# Run on server1:
#   bash scripts/sweep_txcost.sh
#
# Outputs: output/diagnostics/p4_full_txcost/<REGIME>_<tag>.json
#          output/diagnostics/p4_full_txcost/summary.md  (CEV table)

set -e
JULIA=${JULIA:-/home/nanyeon99/.local/bin/julia}
SOLVER=src/vfi_solver_v3.jl
OUTDIR=output/diagnostics/p4_full_txcost
mkdir -p "$OUTDIR"

run_regime () {
    local regime="$1"
    local tag="$2"
    local extra_env="$3"
    local outfile="$OUTDIR/${regime}_${tag}.json"
    if [ -f "$outfile" ]; then
        echo "=== SKIP (exists): regime=$regime tag=$tag ==="
        return
    fi
    echo "=== regime=$regime tag=$tag ==="
    env REGIME="$regime" SUMMARY_JSON_PATH="$outfile" $extra_env \
        "$JULIA" "$SOLVER"
}

for regime in E1_2L E2_2L; do
    run_regime "$regime" "notx"    "TAU_SELL=0.00 TAU_BUY=0.000"
    run_regime "$regime" "sell6"   "TAU_SELL=0.06 TAU_BUY=0.000"
    run_regime "$regime" "rt8p5"   "TAU_SELL=0.06 TAU_BUY=0.025 APPLY_TAU_BUY=1"
    run_regime "$regime" "rt10"    "TAU_SELL=0.06 TAU_BUY=0.040 APPLY_TAU_BUY=1"
    run_regime "$regime" "rt12"    "TAU_SELL=0.06 TAU_BUY=0.060 APPLY_TAU_BUY=1"
done

echo ""
echo "=== Computing CEV sweep table ==="
"$JULIA" scripts/compute_cev_sweep.jl "$OUTDIR" txcost \
    | tee "$OUTDIR/summary.md"
echo ""
echo "Done. Results in $OUTDIR/summary.md"
