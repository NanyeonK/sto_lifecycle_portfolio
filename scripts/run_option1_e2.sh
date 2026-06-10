#!/usr/bin/env bash
# Run E2_2L baseline with Option 1 (6D state, tau_buy on deltas).
# Spec: handoff/tau_buy_option1_spec.md
# Usage: bash scripts/run_option1_e2.sh

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

mkdir -p "$REPO_ROOT/output/diagnostics"

echo "=== Option 1 E2_2L run (vfi_solver_v4.jl) ==="
date

REGIME=E2_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.5 \
ASSET_GRID_SIZE=9 X_GRID_SIZE=5 GH_NODES=3 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
RHO_AB=0.50 \
SUMMARY_JSON_PATH="$REPO_ROOT/output/diagnostics/p6_option1_e2.json" \
JULIA_NUM_THREADS=1 \
  julia "$REPO_ROOT/src/vfi_solver_v4.jl" 2>&1 | tee "$REPO_ROOT/output/diagnostics/p6_option1_e2_stdout.log"

echo "=== E2_2L done. Summary written to output/diagnostics/p6_option1_e2.json ==="
date
