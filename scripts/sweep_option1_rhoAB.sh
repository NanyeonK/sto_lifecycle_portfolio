#!/usr/bin/env bash
# sweep_option1_rhoAB.sh — rho_AB sensitivity sweep for v4 Option-1 solver
# Runs E2_2L at rho_AB ∈ {0, 0.25, 0.50, 0.75, 0.95}; E1_2L baseline unchanged.
# Prerequisite: run_option1_e1.sh must complete first (E1_2L baseline).
# Usage: bash scripts/sweep_option1_rhoAB.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="$REPO_ROOT/output/diagnostics"
mkdir -p "$OUT_DIR"

# Baseline grid settings (match run_option1_e2.sh)
export N_W=15; export N_Z=5; export N_X_PREV=3; export X_PREV_MAX=2.0
export ASSET_GRID_SIZE=7; export X_GRID_SIZE=4; export GH_NODES=3
export SIGMA_DIV=0.10; export GAMMA=5.0
export P_RELOCATE_WORKING=0.06; export P_RELOCATE_RETIRED=0.02
export TAU_SELL=0.06; export TAU_BUY=0.025; export TAU_TOKEN=0.01
export LTV_MAX=0.0
export REGIME=E2_2L

for RHO_AB in 0.00 0.25 0.50 0.75 0.95; do
    SLUG="rhoAB$(echo $RHO_AB | tr -d '.')"
    JSON_OUT="$OUT_DIR/p6_sweep_rhoAB_${SLUG}.json"
    LOG_OUT="$OUT_DIR/p6_sweep_rhoAB_${SLUG}.log"
    echo "=== rho_AB=$RHO_AB starting $(date) ===" | tee "$LOG_OUT"
    RHO_AB=$RHO_AB SUMMARY_JSON_PATH="$JSON_OUT" \
        julia "$REPO_ROOT/src/vfi_solver_v4.jl" 2>&1 | tee -a "$LOG_OUT"
    echo "=== Done $(date) ===" | tee -a "$LOG_OUT"
done

echo "=== rho_AB sweep complete. Files in $OUT_DIR ==="
