#!/usr/bin/env bash
# Run E2_2L baseline under v4 (Option 1 full state extension).
# Produces: output/diagnostics/p6_option1_e2.json
# Run from repo root: bash scripts/run_option1_e2.sh

set -euo pipefail
mkdir -p output/diagnostics

REGIME=E2_2L \
N_W=15 N_Z=5 \
N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 \
GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
RHO=0.05 M=0.01 SIGMA_H=0.115 SIGMA_DIV=0.10 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json \
julia src/vfi_solver_v4.jl

echo "E2_2L v4 done. Summary: output/diagnostics/p6_option1_e2.json"
