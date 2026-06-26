#!/usr/bin/env bash
# Run E1_2L_v4 baseline on server1.
# Uses reduced grids (N_W=15, N_Z=5, N_X_PREV=3) per Option 1 spec.
# Expected wall time: ~2.5 hours single thread.
#
# Usage:  bash scripts/run_option1_e1.sh
# Output: output/diagnostics/p6_option1_e1.json

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

SUMMARY_JSON="$OUTDIR/p6_option1_e1.json"
STDOUT_LOG="$OUTDIR/p6_option1_e1_stdout.log"

echo "=== v4 E1_2L Option 1 baseline ===" | tee "$STDOUT_LOG"
date | tee -a "$STDOUT_LOG"

REGIME=E1_2L \
N_W=15 N_Z=5 \
N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 \
GH_NODES=3 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
SIGMA_S=0.157 G_H=0.016 SIGMA_H=0.115 \
SIGMA_DIV=0.10 \
RHO=0.05 M=0.01 \
SIGMA_U2=0.0106 SIGMA_EPS2=0.0738 LAMBDA_RET=0.65 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
RHO_AB=0.50 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="$SUMMARY_JSON" \
JULIA_NUM_THREADS=1 \
julia src/vfi_solver_v4.jl 2>&1 | tee -a "$STDOUT_LOG"

echo "=== done ===" | tee -a "$STDOUT_LOG"
date | tee -a "$STDOUT_LOG"
echo "Summary JSON: $SUMMARY_JSON"
echo "Stdout log:   $STDOUT_LOG"
