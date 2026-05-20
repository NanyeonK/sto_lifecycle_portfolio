#!/usr/bin/env bash
# Run E2_2L Option 1 (continuous fractional tokens with full x_prev state).
# Output: output/diagnostics/p6_option1_e2.json
# Intended for server1; see handoff/tau_buy_option1_spec.md for context.
#
# Key diagnostic: check mean_xB_t1_xprev00_ellA > 0.
# If mean_xB > 0 at ell=A, x_prev=(0,0): hedge channel is ACTIVE.
# Compare V_t1_midpoint_ellA_xprev00 vs E1_2L run to compute CEV.

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

JULIA_NUM_THREADS=1 \
REGIME=E2_2L \
SUMMARY_JSON_PATH="$OUTDIR/p6_option1_e2.json" \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
SIGMA_S=0.157 SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 \
RHO=0.05 M=0.01 SIGMA_U2=0.0106 SIGMA_EPS2=0.0738 LAMBDA_RET=0.65 \
RHO_AB=0.50 P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 LTV_MAX=0.0 \
AGE0=25 RETIRE_AGE=65 TERMINAL_AGE=80 \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUTDIR/p6_option1_e2_stdout.log"

echo "E2_2L Option 1 run complete. JSON at $OUTDIR/p6_option1_e2.json"
