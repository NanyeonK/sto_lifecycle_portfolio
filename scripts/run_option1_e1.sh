#!/usr/bin/env bash
# run_option1_e1.sh — E1_2L baseline run for v4 (Option 1)
#
# Uses reduced grids to match v3 baseline compute budget:
#   N_W=15, N_Z=5, N_X_PREV=3  → ~4.6x v3 baseline per regime
# Expected wall time: 2-3 hours on server1 (single thread).
#
# Output: output/diagnostics/p6_option1_e1.json + stdout log
#
# Usage: bash scripts/run_option1_e1.sh [--threads N]

set -euo pipefail

THREADS="${1:-1}"
if [[ "$1" == "--threads" && -n "${2:-}" ]]; then THREADS="$2"; fi

mkdir -p output/diagnostics

echo "=== v4 E1_2L baseline (Option 1) ==="
echo "Threads: $THREADS"
echo "Start: $(date)"

JULIA_NUM_THREADS="$THREADS" \
  REGIME=E1_2L \
  N_W=15 \
  N_Z=5 \
  N_X_PREV=3 \
  X_PREV_MAX=1.0 \
  ASSET_GRID_SIZE=7 \
  GH_NODES=3 \
  TAU_SELL=0.06 \
  TAU_BUY=0.025 \
  TAU_TOKEN=0.005 \
  P_RELOCATE_WORKING=0.06 \
  P_RELOCATE_RETIRED=0.02 \
  RHO_AB=0.50 \
  SIGMA_DIV=0.10 \
  GAMMA=5.0 \
  BETA=0.96 \
  SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json \
  julia src/vfi_solver_v4.jl \
  2>&1 | tee output/diagnostics/p6_option1_e1_stdout.log

echo "End: $(date)"
echo "Summary JSON: output/diagnostics/p6_option1_e1.json"
