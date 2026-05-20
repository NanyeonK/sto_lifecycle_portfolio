#!/usr/bin/env bash
# Run E2_2L baseline under Option 1 (v4 solver) on server1.
# E2_2L: continuous fractional tokens of A and/or B; per-period tau_buy on increments.
# The key hedge mechanism: pre-holding x_B at ell=A costs tau_buy incrementally
# but saves full tau_buy on the x_B increment when relocating to B.
#
# Expected wall: ~1-2 hours single thread (6D state; N_W=15, N_Z=5, N_X_PREV=3).
# After run, compare V_t1_midpoint with E1_2L result to compute CEV.
# Key diagnostic: H1_hedge_mean_xB_at_ellA — should be > 0 for mechanism to hold.

set -euo pipefail

mkdir -p output/diagnostics

export REGIME=E2_2L
export N_W=15
export N_Z=5
export N_X_PREV=3
export X_PREV_MAX=1.0
export ASSET_GRID_SIZE=9
export GH_NODES=3

# Baseline calibration (Round 4 confirmed)
export GAMMA=5.0
export BETA=0.96
export RF=1.02
export EQUITY_PREMIUM=0.04
export RHO=0.05
export M=0.01
export SIGMA_H=0.115
export SIGMA_DIV=0.10
export G_H=0.016

# Mobility
export P_RELOCATE_WORKING=0.06
export P_RELOCATE_RETIRED=0.02

# Transaction costs (all active in v4)
export TAU_SELL=0.06
export TAU_BUY=0.025
export TAU_TOKEN=0.01

# Cross-location correlation
export RHO_AB=0.50

export SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json

echo "=== Option 1 E2_2L baseline ==="
echo "  REGIME=$REGIME  N_W=$N_W  N_Z=$N_Z  N_X_PREV=$N_X_PREV"
echo "  tau_sell=$TAU_SELL  tau_buy=$TAU_BUY  tau_token=$TAU_TOKEN"
echo "  Hypothesis 1 test: mean_xB at ell=A should be > 0"
date

julia src/vfi_solver_v4.jl

echo "=== Done ==="
date
