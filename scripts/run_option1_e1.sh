#!/usr/bin/env bash
# Run E1_2L baseline with v4 solver (Option 1: full state extension)
# Execute on server1: bash scripts/run_option1_e1.sh
# Output: output/diagnostics/p6_option1_e1.json

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"

REGIME=E1_2L \
N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 GH_NODES=3 \
TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
SUMMARY_JSON_PATH="$REPO_DIR/output/diagnostics/p6_option1_e1.json" \
julia "$REPO_DIR/src/vfi_solver_v4.jl" 2>&1 | tee "$REPO_DIR/output/diagnostics/p6_option1_e1_stdout.log"

echo "E1_2L v4 done. Summary at output/diagnostics/p6_option1_e1.json"
