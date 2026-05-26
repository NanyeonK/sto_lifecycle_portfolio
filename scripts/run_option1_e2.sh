#!/usr/bin/env bash
# Run E2_2L baseline under v4 solver (Option 1 full state extension)
# Usage: bash scripts/run_option1_e2.sh
# Server1: run AFTER run_option1_e1.sh completes (for CEV comparison)

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "=== v4 Option 1: E2_2L baseline ==="
echo "Grids: N_W=15, N_Z=5, N_X_PREV=3 (small-grid mode)"
echo "Expected wall time: ~2-4 hours (single thread, larger choice space)"

REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=5 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
JULIA_NUM_THREADS=1 \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo "Done. Summary: $OUTDIR/p6_option1_e2.json"
