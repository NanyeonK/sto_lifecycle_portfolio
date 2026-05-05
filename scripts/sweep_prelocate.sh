#!/usr/bin/env bash
# P1 sensitivity: p_relocate_working sweep using v4 solver.
# Sweeps p_relocate_working ∈ {0, 0.02, 0.06, 0.12} for E1_2L and E2_2L.
# At p_relocate=0, the relocation friction disappears and CEV should reflect
# only the continuous-x advantage (no hedge premium).
#
# Run on server1:
#   bash scripts/sweep_prelocate.sh
#
# Outputs: output/diagnostics/p7_prelocate_v4/<REGIME>_preloc<tag>.json
#          output/diagnostics/p7_prelocate_v4/summary.md

set -euo pipefail
JULIA=${JULIA:-julia}
SOLVER=src/vfi_solver_v4.jl
OUTDIR=output/diagnostics/p7_prelocate_v4
mkdir -p "$OUTDIR"

P_RELOC_VALS="0.00 0.02 0.06 0.12"

for p_reloc in $P_RELOC_VALS; do
    tag=$(echo "$p_reloc" | tr '.' 'p')
    for regime in E1_2L E2_2L; do
        outfile="$OUTDIR/${regime}_preloc${tag}.json"
        [ -f "$outfile" ] && { echo "SKIP: $outfile"; continue; }
        echo "=== p_relocate_working=$p_reloc regime=$regime ==="
        REGIME="$regime" \
        P_RELOCATE_WORKING="$p_reloc" \
        N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=2.0 \
        ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
        GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
        SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 RHO=0.05 M=0.01 \
        TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
        P_RELOCATE_RETIRED=0.02 RHO_AB=0.50 LTV_MAX=0.0 \
        SUMMARY_JSON_PATH="$outfile" \
            "$JULIA" "$SOLVER" 2>&1 | tee "${outfile%.json}_stdout.log"
    done
done

echo ""
echo "=== CEV sweep table ==="
"$JULIA" scripts/compute_cev_sweep.jl "$OUTDIR" prelocate | tee "$OUTDIR/summary.md"
echo "Done. $OUTDIR/summary.md"
