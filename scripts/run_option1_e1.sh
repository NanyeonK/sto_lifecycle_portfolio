#!/usr/bin/env bash
# Run E1_2L baseline under Option 1 state extension (v4 solver).
# Output: output/diagnostics/p6_option1_e1.json
# Intended for server1; see handoff/tau_buy_option1_spec.md for context.
#
# Calibration: Round-4 confirmed baseline.
# Grid: N_W=15, N_Z=5, N_X_PREV=3 (coarse; ~2.5h wall estimate).

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

JULIA_NUM_THREADS=1 \
REGIME=E1_2L \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e1.json" \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
SIGMA_S=0.157 SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 \
RHO=0.05 M=0.01 SIGMA_U2=0.0106 SIGMA_EPS2=0.0738 LAMBDA_RET=0.65 \
RHO_AB=0.50 P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 LTV_MAX=0.0 \
AGE0=25 RETIRE_AGE=65 TERMINAL_AGE=80 \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e1_stdout.log"

echo "E1_2L Option 1 run complete. JSON at $OUTDIR/p6_option1_e1.json"
