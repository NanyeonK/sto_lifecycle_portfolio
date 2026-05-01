#!/usr/bin/env bash
# Round 4 P1: p_relocate sensitivity sweep.
# Sweeps p_relocate_working in {0, 0.02, 0.06, 0.12} for E1_2L and E2_2L.
# At p_relocate → 0, cross-location token holding loses its motive (CEV → 0).
#
# Run on server1:
#   bash scripts/sweep_prelocate.sh
#
# Outputs: output/diagnostics/p4_prelocate_sweep/<REGIME>_preloc<tag>.json
#          output/diagnostics/p4_prelocate_sweep/summary.md  (CEV table)

set -e
JULIA=${JULIA:-/home/nanyeon99/.local/bin/julia}
SOLVER=src/vfi_solver_v3.jl
OUTDIR=output/diagnostics/p4_prelocate_sweep
mkdir -p "$OUTDIR"

P_RELOC_VALS="0.00 0.02 0.06 0.12"

for p_reloc in $P_RELOC_VALS; do
    tag=$(echo "$p_reloc" | tr '.' 'p')
    for regime in E1_2L E2_2L; do
        outfile="$OUTDIR/${regime}_preloc${tag}.json"
        if [ -f "$outfile" ]; then
            echo "=== SKIP (exists): p_reloc=$p_reloc regime=$regime ==="
            continue
        fi
        echo "=== p_relocate_working=$p_reloc regime=$regime ==="
        REGIME="$regime" \
        P_RELOCATE_WORKING="$p_reloc" \
        SUMMARY_JSON_PATH="$outfile" \
            "$JULIA" "$SOLVER"
    done
done

echo ""
echo "=== Computing CEV sweep table ==="
"$JULIA" scripts/compute_cev_sweep.jl "$OUTDIR" prelocate \
    | tee "$OUTDIR/summary.md"
echo ""
echo "Done. Results in $OUTDIR/summary.md"
