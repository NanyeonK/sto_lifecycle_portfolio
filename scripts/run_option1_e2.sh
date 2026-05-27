#!/usr/bin/env bash
# run_option1_e2.sh — E2_2L baseline at v4 (6D state, tau_buy active, Option 1)
# Run AFTER run_option1_e1.sh completes and E1_2L JSON is available.
#
# Wall time estimate: ~2.5 hours (single thread, coarse grid)
# Key hypothesis: mean_xB_t1_xp00_ellA > 0  (hedge activates via tau_buy savings)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUT="$REPO_ROOT/output/diagnostics"
mkdir -p "$OUT"

echo "=== Option 1 E2_2L baseline (v4 solver) ==="
echo "repo: $REPO_ROOT"
echo "output: $OUT"
date

REGIME=E2_2L \
  N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
  ASSET_GRID_SIZE=9 GH_NODES=3 \
  TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
  P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
  RHO_AB=0.50 GAMMA=5.0 BETA=0.96 RF=1.02 \
  SIGMA_DIV=0.10 SIGMA_H=0.115 G_H=0.016 \
  LTV_MAX=0.0 \
  SUMMARY_JSON_PATH="$OUT/p6_option1_e2.json" \
  JULIA_NUM_THREADS=1 \
  julia "$REPO_ROOT/src/vfi_solver_v4.jl" \
  2>&1 | tee "$OUT/p6_option1_e2_stdout.log"

echo "=== E2_2L done ===" && date
echo ""
echo "Next: compute CEV(E2_2L vs E1_2L) using scripts/compute_cev_sweep.jl"
echo "  julia scripts/compute_cev_sweep.jl \\"
echo "    output/diagnostics/p6_option1_e1.json \\"
echo "    output/diagnostics/p6_option1_e2.json"
