#!/usr/bin/env bash
# Smoke test for vfi_solver_v4.jl — checks struct allocation, tx_cost, shock block.
# Does NOT run VFI. Safe to run in cloud env (fast, no heavy compute).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

echo "=== Option 1 v4 smoke test ==="
julia src/vfi_solver_v4.jl --smoke-test 2>&1 | tee output/diagnostics/p6_option1_smoke.log

echo ""
echo "Smoke test done. Check output/diagnostics/p6_option1_smoke.log"
