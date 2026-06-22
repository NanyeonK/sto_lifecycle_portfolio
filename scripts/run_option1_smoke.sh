#!/usr/bin/env bash
# Smoke test for v4 solver (fast, no VFI).
# Usage: bash scripts/run_option1_smoke.sh
# Expected wall: < 30 s

set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p output/diagnostics

# Use smallest possible grids for smoke test
export N_X_PREV=3
export X_PREV_MAX=1.0
export GH_NODES=3

julia src/vfi_solver_v4.jl --smoke-test \
    2>&1 | tee output/diagnostics/p6_option1_smoke.md

echo "Smoke test done. See output/diagnostics/p6_option1_smoke.md"
