#!/usr/bin/env bash
# Smoke test v4 solver — no heavy VFI, struct checks + one-period spot-test only.
# Expected wall time: <60 seconds.
set -euo pipefail

mkdir -p output/diagnostics

julia src/vfi_solver_v4.jl --smoke-test 2>&1 | tee output/diagnostics/p6_option1_smoke.log
