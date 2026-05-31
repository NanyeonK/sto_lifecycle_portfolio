#!/usr/bin/env bash
# Falsification tests for Path B Option 1 (v4 solver).
# Under correct model spec, cross-location hedge should ACTIVATE:
#   (r) p_relocate=0 → mean_xB at ellA should = 0 (no relocation, no pre-buy motive)
#   (m) rho_AB=0.95  → hedge should WEAKEN (correlated returns reduce diversification value)
# PASS criteria for Option 1:
#   Baseline (p=0.06, rhoAB=0.5): mean_xB > 0     ← hedge activates
#   p=0.0:                         mean_xB = 0     ← hedge collapses without mobility
#   rhoAB=0.95:                    mean_xB < baseline but possibly > 0 (hedge weakens)
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUT_DIR="$REPO_ROOT/output/diagnostics"
mkdir -p "$OUT_DIR"

BASE_OPTS="N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 ASSET_GRID_SIZE=9 GH_NODES=3 REGIME=E2_2L"
BASE_OPTS="$BASE_OPTS TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 LTV_MAX=0.0"

echo "=== Falsification (r): p_relocate=0.0 ==="
eval "$BASE_OPTS P_RELOCATE_WORKING=0.0 P_RELOCATE_RETIRED=0.0 RHO_AB=0.5 \
  SUMMARY_JSON_PATH=$OUT_DIR/p6_option1_falsify_prelocate0.json \
  julia $REPO_ROOT/src/vfi_solver_v4.jl" \
  2>&1 | tee "$OUT_DIR/p6_option1_falsify_prelocate0_stdout.log"
echo "p_relocate=0 done"

echo "=== Robustness (m): rho_AB=0.95 ==="
eval "$BASE_OPTS P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 RHO_AB=0.95 \
  SUMMARY_JSON_PATH=$OUT_DIR/p6_option1_falsify_rhoAB95.json \
  julia $REPO_ROOT/src/vfi_solver_v4.jl" \
  2>&1 | tee "$OUT_DIR/p6_option1_falsify_rhoAB95_stdout.log"
echo "rho_AB=0.95 done"

echo "=== All falsification tests complete. Check mean_xB_t1_xprev0_ellA in each JSON. ==="
