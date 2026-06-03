#!/usr/bin/env bash
# Run E2_2L baseline with v4 solver (6D state, Option 1 tau_buy).
# Calibration: Round 4 confirmed baseline.
# Expected wall time: ~2–3 h single thread.
# Run on server1: bash scripts/run_option1_e2.sh

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "=== v4 E2_2L baseline (Option 1) ==="
date

REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 \
ASSET_GRID_SIZE=9 GH_NODES=3 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
RHO=0.05 M=0.01 SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 \
RHO_AB=0.50 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo "E2_2L done."
date
