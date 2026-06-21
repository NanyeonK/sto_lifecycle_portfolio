#!/usr/bin/env bash
# Smoke test for v4 solver — no VFI, struct + tx_cost checks only
# Usage: bash scripts/run_option1_smoke.sh
# Expected wall time: < 30 seconds

set -euo pipefail

OUT_DIR="output/diagnostics"
mkdir -p "$OUT_DIR"

JULIA_NUM_THREADS=1 \
julia src/vfi_solver_v4.jl --smoke-test 2>&1 | tee "$OUT_DIR/p6_option1_smoke.log"

echo "Smoke test done. Log at $OUT_DIR/p6_option1_smoke.log"
