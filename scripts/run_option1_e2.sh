#!/usr/bin/env bash
# Run E2_2L with v4 6D-state solver (Option 1 — proper tau_buy hedge via x_prev state).
# Grids: N_W=15, N_Z=5, N_X_PREV=3. Run AFTER run_option1_e1.sh.
# Key hypothesis: mean_xB_t1_xprev0_ellA > 0 (hedge mechanism activates).
# Compare V_t1_midpoint_ellA_xprev0 with E1_2L to compute CEV.
# Usage: bash scripts/run_option1_e2.sh

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=5 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo "E2_2L v4 done. Summary: $OUTDIR/p6_option1_e2.json"
echo ""
echo "Next: compute CEV(E2_2L_v4 vs E1_2L_v4) from the two JSON files."
echo "Hypothesis H1: mean_xB_t1_xprev0_ellA > 0"
echo "Hypothesis H2: CEV > 4.255% (Option 3 baseline)"
