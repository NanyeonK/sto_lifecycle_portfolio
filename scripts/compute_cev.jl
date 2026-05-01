#!/usr/bin/env julia
# compute_cev.jl — compute CEV(E2_2L vs E1_2L) from two JSON summary files.
#
# Usage (from repo root on server1):
#   julia scripts/compute_cev.jl <E1_2L_json> <E2_2L_json> [label]
#
# Outputs: label, V_E1, V_E2, CEV (%) at the representative midpoint state.
#
# CEV formula for CRRA (gamma > 1, negative V):
#   V_E2 = V_E1 * (1 + cev)^(1 - gamma)
#   => cev = (V_E2 / V_E1)^(1 / (1 - gamma)) - 1
# With gamma=5: exponent = 1/(1-5) = -0.25.
# V_E2 > V_E1 (welfare-improving, less negative) => V_E2/V_E1 < 1 => cev > 0. ✓

using JSON3
using Printf

function compute_cev(V1::Float64, V2::Float64, gamma::Float64)::Float64
    V1 >= 0.0 && error("Expected negative V1 (CRRA with gamma>1); got V1=$V1")
    V2 >= 0.0 && error("Expected negative V2; got V2=$V2")
    return (V2 / V1)^(1.0 / (1.0 - gamma)) - 1.0
end

function main(args::Vector{String})
    length(args) < 2 && error("Usage: compute_cev.jl <E1_json> <E2_json> [label]")
    label = length(args) >= 3 ? args[3] : "CEV"

    s1 = JSON3.read(read(args[1], String))
    s2 = JSON3.read(read(args[2], String))

    gamma = Float64(s1["params"]["gamma"])
    V1    = Float64(s1["V_t1_midpoint_ellA"])
    V2    = Float64(s2["V_t1_midpoint_ellA"])

    cev = compute_cev(V1, V2, gamma)

    println("=== $label ===")
    @printf("  regime_E1 : %s\n", s1["regime"])
    @printf("  regime_E2 : %s\n", s2["regime"])
    @printf("  gamma     : %.1f\n", gamma)
    @printf("  V_E1      : %.4f\n", V1)
    @printf("  V_E2      : %.4f\n", V2)
    @printf("  CEV       : %+.4f%%\n", cev * 100.0)

    # Also report mean_xA and mean_xB for the cross-location hedge check
    for (key, s) in [("E1_2L", s1), ("E2_2L", s2)]
        mxA = get(s, "mean_xA_t1_feasible_ellA", nothing)
        mxB = get(s, "mean_xB_t1_feasible_ellA", nothing)
        if mxA !== nothing && mxB !== nothing
            @printf("  %s  mean_xA=%.3f  mean_xB=%.3f\n", key, Float64(mxA), Float64(mxB))
        end
    end
    println()
    return cev
end

main(ARGS)
