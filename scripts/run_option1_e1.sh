#!/usr/bin/env bash
# Run E1_2L baseline under v4 (Option 1 full state extension).
# Usage: bash scripts/run_option1_e1.sh [--full]
# Default: reduced grids (N_W=15, N_Z=5, N_X_PREV=3).
# Pass --full for larger grids (N_W=40, N_Z=9, N_X_PREV=5).

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

if [[ "${1:-}" == "--full" ]]; then
  echo "=== v4 E1_2L FULL GRID ==="
  export N_W=40; export N_Z=9; export N_X_PREV=5; export X_PREV_MAX=2.0
  export ASSET_GRID_SIZE=12; export X_GRID_SIZE=8
  OUT="$OUTDIR/p6_option1_e1_full.json"
  LOG="$OUTDIR/p6_option1_e1_full_stdout.log"
else
  echo "=== v4 E1_2L REDUCED GRID ==="
  export N_W=15; export N_Z=5; export N_X_PREV=3; export X_PREV_MAX=1.5
  export ASSET_GRID_SIZE=7; export X_GRID_SIZE=4
  OUT="$OUTDIR/p6_option1_e1.json"
  LOG="$OUTDIR/p6_option1_e1_stdout.log"
fi

export REGIME=E1_2L
export SUMMARY_JSON_PATH="$OUT"

# Baseline calibration (Round 4 confirmed)
export GAMMA=5.0; export BETA=0.96; export RF=1.02; export EQUITY_PREMIUM=0.04
export SIGMA_H=0.115; export SIGMA_DIV=0.10
export G_H=0.016; export RHO=0.05; export M=0.01
export RHO_AB=0.5
export P_RELOCATE_WORKING=0.06; export P_RELOCATE_RETIRED=0.02
export TAU_SELL=0.06; export TAU_BUY=0.025; export TAU_TOKEN=0.01
export LTV_MAX=0.0; export GH_NODES=3

echo "Output: $OUT"
echo "Log:    $LOG"
time julia src/vfi_solver_v4.jl 2>&1 | tee "$LOG"
echo "Done. Summary at $OUT"
