#!/usr/bin/env bash
# Run E2_2L baseline at v4 settings (Path B Option 1).
# Usage: bash scripts/run_option1_e2.sh
# Run on server1; cloud env may lack Julia.

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"
STAMP=$(date +%Y%m%d_%H%M%S)

echo "=== v4 E2_2L baseline run (Option 1) ==="
echo "  Started: $(date)"

REGIME=E2_2L \
N_W=15 N_Z=5 \
N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 \
GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.5 \
GAMMA=5.0 BETA=0.96 RF=1.02 \
RHO=0.05 M=0.01 SIGMA_H=0.115 SIGMA_DIV=0.10 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2_${STAMP}.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_${STAMP}_stdout.log"

echo "=== E2_2L done. Output: $OUTDIR/p6_option1_e2_${STAMP}.json ==="
