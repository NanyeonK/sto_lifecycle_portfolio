#!/bin/bash
# run_option1_smoke.sh — Quick smoke test for vfi_solver_v4.jl (no VFI run)
# Run first to verify struct-init, tx_cost, and 6D array allocation.
# Should complete in < 5 seconds.
#
# Usage:
#   bash scripts/run_option1_smoke.sh

set -euo pipefail

OUTPUT_DIR="output/diagnostics"
mkdir -p "$OUTPUT_DIR"

LOG="$OUTPUT_DIR/p6_option1_smoke.md"

echo "=== v4 smoke test: $(date) ===" | tee "$LOG"

julia src/vfi_solver_v4.jl --smoke-test 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== smoke test done: $(date) ===" | tee -a "$LOG"
echo "Log at: $LOG"
