#!/usr/bin/env bash
# Smoke test for v4 solver — struct checks only, no VFI.
# Safe to run in cloud or on server1.

set -euo pipefail

julia src/vfi_solver_v4.jl --smoke-test
