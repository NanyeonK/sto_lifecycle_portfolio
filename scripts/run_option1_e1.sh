#!/usr/bin/env bash
# Run E1_2L baseline under v4 (Option 1 full state extension).
# Grid: N_W=15, N_Z=5, N_X_PREV=3 (default small mode).
# Wall estimate: ~1-2 hours single thread.
set -euo pipefail

REGIME=E1_2L
OUT=output/diagnostics/p6_option1_e1
mkdir -p output/diagnostics

echo "=== v4 E1_2L baseline (Option 1) ==="
REGIME=${REGIME} \
SUMMARY_JSON_PATH=${OUT}.json \
julia src/vfi_solver_v4.jl 2>&1 | tee ${OUT}_stdout.log

echo "Done. Summary at ${OUT}.json"
echo "stdout log at ${OUT}_stdout.log"
