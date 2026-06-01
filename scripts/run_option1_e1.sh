#!/usr/bin/env bash
# Run E1_2L baseline at v4 (Option 1) grid settings.
# Expected wall time: ~2-3 hours single thread on server1.
# Output: output/diagnostics/p6_option1_e1.json

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
mkdir -p output/diagnostics

REGIME=E1_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=7 X_GRID_SIZE=4 \
GH_NODES=3 \
GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04 \
SIGMA_S=0.157 G_H=0.016 SIGMA_H=0.115 SIGMA_DIV=0.10 \
RHO=0.05 M=0.01 \
RHO_AB=0.5 P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json \
  julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e1_stdout.log

echo "Done. Summary: output/diagnostics/p6_option1_e1.json"
