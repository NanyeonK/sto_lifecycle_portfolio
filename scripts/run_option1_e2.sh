#!/usr/bin/env bash
# run_option1_e2.sh — E2_2L baseline run for v4 (Option 1)
#
# Same grid as E1_2L for direct CEV comparison.
# E2_2L has N_X_PREV^2 = 9 housing choices per state (vs 2 for E1_2L),
# so actual per-period compute is ~4.5x E1_2L; budget accordingly.
#
# Expected wall time: ~10-12 hours on server1 (single thread).
# With --threads 4: ~3 hours.
#
# Output: output/diagnostics/p6_option1_e2.json + stdout log
#
# Usage: bash scripts/run_option1_e2.sh [--threads N]

set -euo pipefail

THREADS=1
if [[ "${1:-}" == "--threads" && -n "${2:-}" ]]; then THREADS="$2"; fi

mkdir -p output/diagnostics

echo "=== v4 E2_2L baseline (Option 1) ==="
echo "Threads: $THREADS"
echo "Start: $(date)"

JULIA_NUM_THREADS="$THREADS" \
  REGIME=E2_2L \
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
  SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json \
  julia src/vfi_solver_v4.jl \
  2>&1 | tee output/diagnostics/p6_option1_e2_stdout.log

echo "End: $(date)"
echo "Summary JSON: output/diagnostics/p6_option1_e2.json"
