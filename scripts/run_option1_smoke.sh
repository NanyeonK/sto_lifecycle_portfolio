#!/usr/bin/env bash
# Smoke test for v4 solver (no VFI; runs in seconds).
# Execute on server1: bash scripts/run_option1_smoke.sh
set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "=== v4 smoke test ==="
date

julia src/vfi_solver_v4.jl --smoke-test 2>&1 | tee "$OUTDIR/p6_option1_smoke_stdout.log"

echo "Smoke test done."
date
