#!/usr/bin/env bash
# Run E2_2L baseline at v4 (Option 1) settings.
# Key test: does mean_xB > 0 at ell=A?  (hedge channel activation)
# Expected wall time: ~2.5h single thread.
# Run on server1 in tmux session sto_lifecycle_portfolio.

set -euo pipefail

OUT="output/diagnostics"
mkdir -p "$OUT"

SUMMARY_JSON_PATH="$OUT/p6_option1_e2.json" \
REGIME=E2_2L \
N_W=15  N_Z=5  N_X_PREV=3  X_PREV_MAX=2.0 \
ASSET_GRID_SIZE=9  X_GRID_SIZE=5  GH_NODES=3 \
GAMMA=5.0  BETA=0.96  RF=1.02  EQUITY_PREMIUM=0.04 \
RHO=0.05   M=0.01     SIGMA_H=0.115  SIGMA_DIV=0.10 \
G_H=0.016  RHO_AB=0.50 \
P_RELOCATE_WORKING=0.06  P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06  TAU_BUY=0.025  TAU_TOKEN=0.01 \
LTV_MAX=0.0 \
julia src/vfi_solver_v4.jl 2>&1 | tee "$OUT/p6_option1_e2_stdout.log"

echo "E2_2L v4 done. Summary: $OUT/p6_option1_e2.json"
