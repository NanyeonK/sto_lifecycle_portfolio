#!/usr/bin/env bash
# Run E1_2L baseline with v4 solver (Option 1 full state extension)
# Usage: bash scripts/run_option1_e1.sh
# Expected wall time: ~2-3 hours single thread (6D state, N_X_PREV=3)
#
# X_PREV_MAX=1.0 is required so that x=1.0 (binary "own") falls exactly on the
# x_prev grid {0, 0.5, 1.0}. Using X_PREV_MAX=1.5 would put the own state (x=1)
# off-grid (nearest neighbor 0.75 ≠ 1) and break E1_2L binary ownership logic.

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

JULIA="${JULIA_PATH:-julia}"

echo "=== v4 E1_2L baseline (Option 1) ==="
echo "Start: $(date)"

REGIME=E1_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
RHO=0.05 M=0.01 SIGMA_H=0.115 SIGMA_DIV=0.10 LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e1.json" \
JULIA_NUM_THREADS=1 \
"$JULIA" src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e1_stdout.log"

echo "Done: $(date)"
