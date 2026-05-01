#!/usr/bin/env julia
# compute_cev_v3.jl — compute CEV across v3 solver JSON summary pairs
#
# Usage:
#   julia scripts/compute_cev_v3.jl <baseline_json> <alternative_json>
#
# Example (Round 4 channel decomposition):
#   julia scripts/compute_cev_v3.jl \
#       output/diagnostics/p4_txcost_E1_2L_baseline.json \
#       output/diagnostics/p4_txcost_E2_2L_baseline.json
#
# CEV formula (CRRA, gamma ≠ 1):
#   V_alt = V_base * (1 + CEV)^(1-gamma)  at each state
#   CEV = (V_alt / V_base)^(1/(1-gamma)) - 1
#
# Reports CEV at midpoint state and mean over all feasible t=1 states.

using JSON3

function read_summary(path::String)
    open(path) do io; JSON3.read(io, Dict{String,Any}); end
end

function compute_cev(V_base::Float64, V_alt::Float64, gamma::Float64)
    (V_base >= 0.0 || V_alt >= 0.0) && return NaN   # non-negative values unexpected
    !isfinite(V_base) || !isfinite(V_alt) && return NaN
    ratio = V_alt / V_base
    ratio <= 0.0 && return NaN
    return ratio^(1.0 / (1.0 - gamma)) - 1.0
end

function main()
    length(ARGS) < 2 && error("Usage: compute_cev_v3.jl <base.json> <alt.json>")
    base_path, alt_path = ARGS[1], ARGS[2]

    base = read_summary(base_path)
    alt  = read_summary(alt_path)

    gamma = Float64(base["params"]["gamma"])
    @assert Float64(alt["params"]["gamma"]) == gamma "gamma mismatch between scenarios"

    println("CEV computation: $(base["regime"]) → $(alt["regime"])")
    println("  gamma = $gamma")
    println("  base file: $base_path")
    println("  alt  file: $alt_path")
    println()

    V_base_mid = Float64(base["V_t1_midpoint_ellA"])
    V_alt_mid  = Float64(alt["V_t1_midpoint_ellA"])
    cev_mid    = compute_cev(V_base_mid, V_alt_mid, gamma)
    @printf("  CEV at midpoint state (ellA): %+.3f%%\n", 100.0 * cev_mid)

    V_base_mean = Float64(base["V_t1_mean_feasible_ellA"])
    V_alt_mean  = Float64(alt["V_t1_mean_feasible_ellA"])
    cev_mean    = compute_cev(V_base_mean, V_alt_mean, gamma)
    @printf("  CEV at mean feasible state (ellA): %+.3f%%\n", 100.0 * cev_mean)

    println()
    println("  [ellB symmetric check]")
    V_base_midB = Float64(base["V_t1_midpoint_ellB"])
    V_alt_midB  = Float64(alt["V_t1_midpoint_ellB"])
    cev_midB    = compute_cev(V_base_midB, V_alt_midB, gamma)
    @printf("  CEV at midpoint state (ellB): %+.3f%%\n", 100.0 * cev_midB)
end

main()
