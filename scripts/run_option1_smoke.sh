#!/usr/bin/env bash
# Smoke test for vfi_solver_v4.jl — no VFI run, checks struct + logic only
# Run this first to verify the v4 implementation before running full VFI.
# Expected wall time: < 30 seconds (no quadrature or VFI).
#
# Usage:
#   bash scripts/run_option1_smoke.sh

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$REPO/output/diagnostics"
mkdir -p "$OUT"

echo "=== Option 1 v4 smoke test ==="
echo "Branch: $(git -C "$REPO" rev-parse --abbrev-ref HEAD)"
echo "Commit: $(git -C "$REPO" rev-parse --short HEAD)"
echo ""

julia "$REPO/src/vfi_solver_v4.jl" --smoke-test 2>&1 | tee "$OUT/p6_option1_smoke.log"

echo ""
echo "Smoke log written to: $OUT/p6_option1_smoke.log"
