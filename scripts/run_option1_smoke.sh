#!/usr/bin/env bash
# Smoke test for v4 solver (no VFI; structural checks only).
# Run from repo root: bash scripts/run_option1_smoke.sh
set -euo pipefail
julia src/vfi_solver_v4.jl --smoke-test
