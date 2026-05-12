#!/usr/bin/env bash
# Run E1_2L baseline under v4 Option 1 (6D state, per-period tx_cost).
# Usage: bash scripts/run_option1_e1.sh
# Output written to output/diagnostics/p6_option1_e1.json

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "=== v4 E1_2L baseline ==="
REGIME=E1_2L \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
N_W=15 \
N_Z=5 \
ASSET_GRID_SIZE=7 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5.0 \
BETA=0.96 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e1.json" \
julia src/vfi_solver_v4.jl

echo "Done. Summary at $OUTDIR/p6_option1_e1.json"
