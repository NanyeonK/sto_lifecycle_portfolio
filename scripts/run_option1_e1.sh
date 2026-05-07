#!/usr/bin/env bash
# Run E1_2L baseline under Option 1 (v4 solver with x_prev state).
# Matches Option 1 spec calibration. Wall time: ~30-90 min on server1.
#
# Usage: bash scripts/run_option1_e1.sh
# Prerequisite: Julia available at julia or $JULIA_BIN; JSON3 installed.

set -euo pipefail

JULIA_BIN="${JULIA_BIN:-julia}"
OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "=== Option 1 — E1_2L baseline (v4 solver) ==="
echo "Started: $(date)"

REGIME=E1_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 GH_NODES=3 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
RHO=0.05 M=0.01 SIGMA_H=0.115 SIGMA_DIV=0.10 \
RHO_AB=0.50 LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e1.json" \
"$JULIA_BIN" src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e1_stdout.log"

echo "Finished: $(date)"
echo "Summary: $OUTDIR/p6_option1_e1.json"
