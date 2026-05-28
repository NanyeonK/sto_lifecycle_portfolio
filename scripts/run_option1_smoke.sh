#!/usr/bin/env bash
# Smoke test for v4 solver — structural checks only, no VFI.
# Run on server1: bash scripts/run_option1_smoke.sh
set -euo pipefail

echo "=== v4 smoke test ==="
julia src/vfi_solver_v4.jl --smoke-test 2>&1 | tee output/diagnostics/p6_option1_smoke.log
echo "Smoke test complete. Results in output/diagnostics/p6_option1_smoke.log"
