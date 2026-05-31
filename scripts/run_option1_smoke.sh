#!/usr/bin/env bash
# Smoke test for v4 solver (Option 1: 6D state extension).
# Run on server1 to verify struct layout, tx_cost logic, and shock block.
# No VFI is executed — safe to run in any environment.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
julia "$REPO_ROOT/src/vfi_solver_v4.jl" --smoke-test
