#!/usr/bin/env bash
# Round 4 P1: rho_AB sensitivity sweep.
# Sweeps rho_AB in {0, 0.25, 0.50, 0.75, 0.95} for both E1_2L and E2_2L.
# At rho_AB → 1, the cross-location hedge channel should collapse (CEV → 0).
#
# Run on server1:
#   bash scripts/sweep_rhoAB.sh
#
# Outputs: output/diagnostics/p4_rhoAB_sweep/<REGIME>_rhoAB<tag>.json
#          output/diagnostics/p4_rhoAB_sweep/summary.md  (CEV table)

set -e
JULIA=${JULIA:-/home/nanyeon99/.local/bin/julia}
SOLVER=src/vfi_solver_v3.jl
OUTDIR=output/diagnostics/p4_rhoAB_sweep
mkdir -p "$OUTDIR"

RHO_AB_VALS="0.00 0.25 0.50 0.75 0.95"

for rho_ab in $RHO_AB_VALS; do
    tag=$(echo "$rho_ab" | tr '.' 'p')
    for regime in E1_2L E2_2L; do
        outfile="$OUTDIR/${regime}_rhoAB${tag}.json"
        if [ -f "$outfile" ]; then
            echo "=== SKIP (exists): rho_AB=$rho_ab regime=$regime ==="
            continue
        fi
        echo "=== rho_AB=$rho_ab regime=$regime ==="
        REGIME="$regime" \
        RHO_AB="$rho_ab" \
        SUMMARY_JSON_PATH="$outfile" \
            "$JULIA" "$SOLVER"
    done
done

echo ""
echo "=== Computing CEV sweep table ==="
"$JULIA" scripts/compute_cev_sweep.jl "$OUTDIR" rhoAB \
    | tee "$OUTDIR/summary.md"
echo ""
echo "Done. Results in $OUTDIR/summary.md"
