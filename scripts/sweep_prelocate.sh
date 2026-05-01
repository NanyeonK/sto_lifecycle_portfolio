#!/usr/bin/env bash
# sweep_prelocate.sh — Round 4 P1 sensitivity: p_relocate_working in {0, 0.02, 0.06, 0.12}
#
# Runs E1_2L and E2_2L at each mobility rate. Writes per-run JSON summaries.
# Final output: output/diagnostics/p4_prelocate_sweep.md assembled by this script.
#
# Usage: bash scripts/sweep_prelocate.sh
# Requires: Julia at $JULIA (default: julia), vfi_solver_v3.jl in src/
#
# Expected mechanism prediction:
#   At p_relocate = 0: no moves -> relocation friction irrelevant ->
#     E2_2L tokens provide no advantage over E1_2L -> CEV -> 0
#     (cross-location holding collapses; no relocation events)
#   At p_relocate = 0.12: high mobility -> relocation frequent ->
#     transaction costs hit often in E1_2L -> E2_2L advantage largest
#   PSID working-age baseline: 0.06 (default)

set -euo pipefail

JULIA="${JULIA:-julia}"
SOLVER="src/vfi_solver_v3.jl"
OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

P_RELOCATE_VALUES="0 0.02 0.06 0.12"

SUMMARY_FILE="$OUTDIR/p4_prelocate_sweep.md"
echo "# Round 4 P1: p_relocate Sensitivity Sweep" > "$SUMMARY_FILE"
echo "Date: $(date -u +%Y-%m-%d)" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "| p_relocate | Regime | V_t1_midpoint_ellA | mean_xA_ellA | mean_xB_ellA | xB_gt0_count |" >> "$SUMMARY_FILE"
echo "|------------|--------|--------------------|--------------|--------------|--------------|" >> "$SUMMARY_FILE"

for p_reloc in $P_RELOCATE_VALUES; do
    p_tag=$(echo "$p_reloc" | tr '.' 'p')
    for regime in E1_2L E2_2L; do
        json_out="$OUTDIR/p4_prelocate_${p_tag}_${regime}.json"
        log_out="$OUTDIR/p4_prelocate_${p_tag}_${regime}_stdout.log"

        echo "--- p_relocate_working=$p_reloc  regime=$regime ---"
        REGIME="$regime" \
        P_RELOCATE_WORKING="$p_reloc" \
        SUMMARY_JSON_PATH="$json_out" \
        "$JULIA" --threads=1 "$SOLVER" 2>&1 | tee "$log_out"

        if [ -f "$json_out" ]; then
            V=$(python3 -c "import json,sys; d=json.load(open('$json_out')); print(round(d.get('V_t1_midpoint_ellA', float('nan')), 4))" 2>/dev/null || echo "N/A")
            mxA=$(python3 -c "import json,sys; d=json.load(open('$json_out')); print(round(d.get('mean_xA_t1_feasible_ellA', float('nan')), 4))" 2>/dev/null || echo "N/A")
            mxB=$(python3 -c "import json,sys; d=json.load(open('$json_out')); print(round(d.get('mean_xB_t1_feasible_ellA', float('nan')), 4))" 2>/dev/null || echo "N/A")
            cnt=$(python3 -c "import json,sys; d=json.load(open('$json_out')); print(d.get('xB_gt0_count_t1_ellA', 'N/A'))" 2>/dev/null || echo "N/A")
            echo "| $p_reloc | $regime | $V | $mxA | $mxB | $cnt |" >> "$SUMMARY_FILE"
        else
            echo "| $p_reloc | $regime | (run failed) | | | |" >> "$SUMMARY_FILE"
        fi
    done
done

# Append interpretation note
cat >> "$SUMMARY_FILE" << 'INTERP'

## Interpretation

Mobility prediction: `CEV(E2_2L vs E1_2L)` should be monotone increasing in `p_relocate`.
At `p_relocate = 0`: no relocation -> no transaction cost friction in E1_2L -> no E2_2L advantage.
Cross-location xB holdings should also collapse at p_relocate=0 because there is no relocation risk to hedge.
At `p_relocate = 0.12`: double the PSID baseline -> more frequent forced sales -> larger E2_2L benefit.

PSID anchor: working-age household-level mobility rate ~5-7% per year. Baseline: 0.06.

CEV computation: use the paired V values from this table to compute
`CEV(E2_2L vs E1_2L) = (V_E2 / V_E1)^(1/(1-gamma)) - 1` for each p_relocate.
INTERP

echo ""
echo "=== p_relocate sweep complete. Summary: $SUMMARY_FILE ==="
