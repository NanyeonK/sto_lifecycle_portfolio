#!/usr/bin/env bash
# Run E1_2L baseline with v4 Option 1 state extension.
# Uses coarse x_prev grid (N_X_PREV=3, X_PREV_MAX=1.0) and reduced N_W, N_Z.
# Expected wall time: ~2.5h single thread at default grid settings.
# Run from repo root. Julia must be in PATH.
set -euo pipefail

mkdir -p output/diagnostics

REGIME=E1_2L \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
N_W=15 \
N_Z=5 \
ASSET_GRID_SIZE=9 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5 \
SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json \
JULIA_NUM_THREADS=1 \
julia src/vfi_solver_v4.jl 2>&1 | tee output/diagnostics/p6_option1_e1_stdout.log

echo "E1_2L v4 done. Summary: output/diagnostics/p6_option1_e1.json"
