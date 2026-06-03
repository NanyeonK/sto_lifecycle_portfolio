#!/usr/bin/env bash
# Smoke test for v4 solver — no VFI, struct/shock/array checks only.
# Fast to run; safe for cloud env.
set -euo pipefail

julia src/vfi_solver_v4.jl --smoke-test 2>&1 | tee output/diagnostics/p6_option1_smoke.log
echo "Smoke test done."
