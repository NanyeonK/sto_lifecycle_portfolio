#!/usr/bin/env bash
# sweep_rhoAB.sh — Round 4 P1 sensitivity: rho_AB in {0, 0.25, 0.5, 0.75, 0.95}
#
# Runs E1_2L and E2_2L at each rho_AB value. Writes per-run JSON summaries.
# Final output: output/diagnostics/p4_rhoAB_sweep.md assembled by this script.
#
# Usage: bash scripts/sweep_rhoAB.sh
# Requires: Julia at $JULIA (default: julia), vfi_solver_v3.jl in src/
#
# Expected hedge-channel prediction:
#   At rho_AB -> 1: locations perfectly correlated -> E2_2L hedge value -> 0
#   At rho_AB  = 0: orthogonal markets -> maximum hedge benefit
#   Maintained-hedge channel = CEV(E2_2L vs E1_2L) should decline with rho_AB.

set -euo pipefail

JULIA="${JULIA:-julia}"
SOLVER="src/vfi_solver_v3.jl"
OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

RHO_AB_VALUES="0 0.25 0.5 0.75 0.95"

SUMMARY_FILE="$OUTDIR/p4_rhoAB_sweep.md"
echo "# Round 4 P1: rho_AB Sensitivity Sweep" > "$SUMMARY_FILE"
echo "Date: $(date -u +%Y-%m-%d)" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "| rho_AB | Regime | V_t1_midpoint_ellA | mean_xA_ellA | mean_xB_ellA | xB_gt0_count |" >> "$SUMMARY_FILE"
echo "|--------|--------|--------------------|--------------|--------------|--------------|" >> "$SUMMARY_FILE"

for rho_ab in $RHO_AB_VALUES; do
    rho_tag=$(echo "$rho_ab" | tr '.' 'p')
    for regime in E1_2L E2_2L; do
        json_out="$OUTDIR/p4_rhoAB_${rho_tag}_${regime}.json"
        log_out="$OUTDIR/p4_rhoAB_${rho_tag}_${regime}_stdout.log"

        echo "--- rho_AB=$rho_ab  regime=$regime ---"
        REGIME="$regime" \
        RHO_AB="$rho_ab" \
        SUMMARY_JSON_PATH="$json_out" \
        "$JULIA" --threads=1 "$SOLVER" 2>&1 | tee "$log_out"

        # Extract key fields from JSON and append table row
        if [ -f "$json_out" ]; then
            V=$(python3 -c "import json,sys; d=json.load(open('$json_out')); print(round(d.get('V_t1_midpoint_ellA', float('nan')), 4))" 2>/dev/null || echo "N/A")
            mxA=$(python3 -c "import json,sys; d=json.load(open('$json_out')); print(round(d.get('mean_xA_t1_feasible_ellA', float('nan')), 4))" 2>/dev/null || echo "N/A")
            mxB=$(python3 -c "import json,sys; d=json.load(open('$json_out')); print(round(d.get('mean_xB_t1_feasible_ellA', float('nan')), 4))" 2>/dev/null || echo "N/A")
            cnt=$(python3 -c "import json,sys; d=json.load(open('$json_out')); print(d.get('xB_gt0_count_t1_ellA', 'N/A'))" 2>/dev/null || echo "N/A")
            echo "| $rho_ab | $regime | $V | $mxA | $mxB | $cnt |" >> "$SUMMARY_FILE"
        else
            echo "| $rho_ab | $regime | (run failed) | | | |" >> "$SUMMARY_FILE"
        fi
    done
done

# Append interpretation note
cat >> "$SUMMARY_FILE" << 'INTERP'

## Interpretation

Hedge-channel prediction: `CEV(E2_2L vs E1_2L)` should decline monotonically with `rho_AB`.
At `rho_AB -> 1`, E2_2L cross-location holding loses value because A and B move together —
the household cannot diversify by retaining the other location's token.
At `rho_AB = 0`, the two markets are orthogonal and the hedge channel is maximized.

CEV computation: use the paired V values from this table to compute
`CEV(E2_2L vs E1_2L) = (V_E2 / V_E1)^(1/(1-gamma)) - 1` for each rho_AB.
(Requires same grid and gamma for comparability.)
INTERP

echo ""
echo "=== rho_AB sweep complete. Summary: $SUMMARY_FILE ==="
