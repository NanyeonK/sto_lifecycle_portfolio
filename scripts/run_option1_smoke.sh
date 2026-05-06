#!/usr/bin/env bash
# Smoke test for v4 solver — struct init and tx_cost checks only; no VFI.
# Run on server1 first before the full baseline runs.
# Expected wall time: < 30 seconds.

set -euo pipefail
cd "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"

OUTDIR=output/diagnostics
mkdir -p "$OUTDIR"

julia src/vfi_solver_v4.jl --smoke-test 2>&1 | tee "$OUTDIR/p6_option1_smoke.md"

echo "Smoke test complete. Output: $OUTDIR/p6_option1_smoke.md"
