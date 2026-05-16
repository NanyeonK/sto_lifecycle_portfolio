#!/usr/bin/env bash
# sweep_option1_prelocate.sh — p_relocate sensitivity sweep for v4 Option-1 solver
# Runs BOTH E1_2L and E2_2L at p_relocate_working ∈ {0, 0.06, 0.12, 0.30}.
# At p_relocate=0: hedge channel must collapse (mean_xB → 0 in E2_2L).
# At p_relocate=0.30: high mobility, hedge premium maximized.
# Usage: bash scripts/sweep_option1_prelocate.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="$REPO_ROOT/output/diagnostics"
mkdir -p "$OUT_DIR"

export N_W=15; export N_Z=5; export N_X_PREV=3; export X_PREV_MAX=2.0
export ASSET_GRID_SIZE=7; export X_GRID_SIZE=4; export GH_NODES=3
export SIGMA_DIV=0.10; export RHO_AB=0.50; export GAMMA=5.0
export TAU_SELL=0.06; export TAU_BUY=0.025; export TAU_TOKEN=0.01
export LTV_MAX=0.0

for P_RELOC in 0.00 0.06 0.12 0.30; do
    SLUG="preloc$(echo $P_RELOC | tr -d '.')"
    for REGIME in E1_2L E2_2L; do
        JSON_OUT="$OUT_DIR/p6_sweep_prelocate_${REGIME}_${SLUG}.json"
        LOG_OUT="$OUT_DIR/p6_sweep_prelocate_${REGIME}_${SLUG}.log"
        echo "=== REGIME=$REGIME p_reloc=$P_RELOC starting $(date) ===" | tee "$LOG_OUT"
        REGIME=$REGIME P_RELOCATE_WORKING=$P_RELOC P_RELOCATE_RETIRED=$(echo "$P_RELOC * 0.33" | bc -l) \
            SUMMARY_JSON_PATH="$JSON_OUT" \
            julia "$REPO_ROOT/src/vfi_solver_v4.jl" 2>&1 | tee -a "$LOG_OUT"
        echo "=== Done $(date) ===" | tee -a "$LOG_OUT"
    done
done

echo "=== p_relocate sweep complete. Files in $OUT_DIR ==="
