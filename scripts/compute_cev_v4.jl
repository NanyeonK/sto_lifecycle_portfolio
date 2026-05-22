#!/usr/bin/env julia
# compute_cev_v4.jl — CEV computation for v4 Option 1 results.
#
# Usage (single pair, H1/H2/H3 verdict):
#   julia scripts/compute_cev_v4.jl baseline \
#         output/diagnostics/p6_option1_e1.json \
#         output/diagnostics/p6_option1_e2.json
#
# Usage (sweep table, after running sweep_v4_rhoAB.sh or sweep_v4_prelocate.sh):
#   julia scripts/compute_cev_v4.jl rhoAB    output/diagnostics/p6_v4_rhoAB/
#   julia scripts/compute_cev_v4.jl prelocate output/diagnostics/p6_v4_prelocate/
#   julia scripts/compute_cev_v4.jl mortgage  output/diagnostics/p6_v4_mortgage/
#
# CEV formula: C = (V_E2 / V_E1)^(1/(1-γ)) - 1.
# For γ=5: C = (V_E2/V_E1)^(-1/4) - 1. Both V < 0; V_E2 > V_E1 → ratio < 1 → C > 0.
#
# v4 key differences from v3 summary JSON:
#   - Representative V key: "V_t1_midpoint_ellA_prevzero" (x_prev=(0,0) slice)
#   - New hedge-activation key: "xB_positive_frac_t1_ellA"
#   - No "apply_tau_buy_at_reloc" key.

using JSON3, Printf, Dates

function compute_cev(V_E2::Float64, V_E1::Float64, gamma::Float64)::Float64
    ratio = V_E2 / V_E1
    ratio <= 0.0 && error("V_E2/V_E1 = $ratio ≤ 0 — check sign of value function")
    return ratio^(1.0 / (1.0 - gamma)) - 1.0
end

function read_v4_summary(path::String)
    isfile(path) || error("File not found: $path")
    d = JSON3.read(read(path, String))
    V_rep     = Float64(d["V_t1_midpoint_ellA_prevzero"])
    gm        = Float64(d["params"]["gamma"])
    rhoAB     = Float64(d["params"]["rho_AB"])
    preloc    = Float64(d["params"]["p_relocate_working"])
    ts        = Float64(d["params"]["tau_sell"])
    tb        = Float64(d["params"]["tau_buy"])
    tt        = Float64(d["params"]["tau_token"])
    ltv       = Float64(d["params"]["ltv_max"])
    # Hedge-channel indicators
    mean_xB   = haskey(d, "mean_xB_t1_feasible_ellA") && d["mean_xB_t1_feasible_ellA"] !== nothing ?
                Float64(d["mean_xB_t1_feasible_ellA"]) : NaN
    xB_frac   = haskey(d, "xB_positive_frac_t1_ellA") && d["xB_positive_frac_t1_ellA"] !== nothing ?
                Float64(d["xB_positive_frac_t1_ellA"]) : NaN
    return (; V_rep, gm, rhoAB, preloc, ts, tb, tt, ltv, mean_xB, xB_frac, path)
end

function print_baseline_verdict(path_e1::String, path_e2::String)
    println("# CEV Verdict: Option 1 v4 Baseline")
    println("Generated: $(Dates.now())")
    println()

    s1 = read_v4_summary(path_e1)
    s2 = read_v4_summary(path_e2)
    cev = compute_cev(s2.V_rep, s1.V_rep, s1.gm)

    println("## Representative state (midpoint w, z | ell=A | x_prev=(0,0))")
    @printf("  V(E1_2L): %.4f\n", s1.V_rep)
    @printf("  V(E2_2L): %.4f\n", s2.V_rep)
    @printf("  CEV(E2_2L vs E1_2L_v4): %+.3f%%\n", cev * 100)
    println()

    println("## H1/H2/H3 Hypothesis Verdicts")
    # H1: mean_xB > 0 at ellA (hedge mechanism activates)
    h1_pass = !isnan(s2.mean_xB) && s2.mean_xB > 1e-4
    @printf("  H1: mean_xB at ellA = %.4f  → %s\n",
            isnan(s2.mean_xB) ? 0.0 : s2.mean_xB,
            h1_pass ? "PASS (hedge activates)" : "FAIL (hedge dead)")
    if !isnan(s2.xB_frac)
        @printf("       xB_positive_frac = %.3f%%\n", s2.xB_frac * 100)
    end

    # H2: CEV > 4.255% (Option 3 baseline)
    option3_baseline = 0.04255
    h2_pass = cev > option3_baseline
    @printf("  H2: CEV = %+.3f%%  vs Option 3 baseline 4.255%%  → %s\n",
            cev * 100,
            h2_pass ? "PASS (above baseline)" : "FAIL (below baseline)")

    # H3: incremental hedge channel (require both E2_2L_v3 baseline ~3.995% for comparison)
    println("  H3: CEV(E2_2L_v4 vs E2_2L_v3) ≈ 0.5-1.5%  → compute manually")
    println("      Reference: V(E2_2L_v3) at same calibration (~-1193.5 from prior logs)")
    println()

    println("## Calibration Used")
    @printf("  gamma=%.1f, rho_AB=%.2f, p_relocate_working=%.2f\n",
            s1.gm, s1.rhoAB, s1.preloc)
    @printf("  tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f, ltv_max=%.2f\n",
            s1.ts, s1.tb, s1.tt, s1.ltv)
    println()

    all_pass = h1_pass && h2_pass
    println("## Decision")
    if all_pass
        println("  H1+H2: PASS → v4 mechanism credible. Run P1 sensitivity sweeps.")
        println("  If H3 also in 0.5-1.5% range: RFS-marginal. Continue to Phase 2.")
    elseif h2_pass && !h1_pass
        println("  H2: PASS but H1: FAIL — CEV gain exists but not from hedge.")
        println("  Mechanism is still continuous-x + tx-cost (same as v3 at higher CEV).")
        println("  Decision: verify channel decomposition before RFS claim.")
    else
        println("  H1+H2 FAIL → v4 does not improve over Option 3 baseline.")
        println("  Fall back to PATH D: REE/JHE with current +4.26% evidence.")
    end
end

function cev_row_v4(dir, tag_e1, tag_e2)
    f1 = joinpath(dir, tag_e1); f2 = joinpath(dir, tag_e2)
    (isfile(f1) && isfile(f2)) || return nothing, "MISSING: $(basename(f1)) or $(basename(f2))"
    s1 = read_v4_summary(f1); s2 = read_v4_summary(f2)
    cev = compute_cev(s2.V_rep, s1.V_rep, s1.gm)
    return (; s1, s2, cev), nothing
end

function print_sweep_rhoAB(dir)
    println("# rho_AB Sensitivity — v4 Option 1")
    println("Generated: $(Dates.now())")
    println()
    println("| rho_AB | V(E1_2L) | V(E2_2L) | CEV (%) | mean_xB@A | xB_frac | note |")
    println("|--------|----------|----------|---------|-----------|---------|------|")
    for (val, tag) in [(0.00,"0p00"), (0.25,"0p25"), (0.50,"0p50"), (0.75,"0p75"), (0.95,"0p95")]
        row, err = cev_row_v4(dir, "E1_2L_rhoAB$(tag).json", "E2_2L_rhoAB$(tag).json")
        if row === nothing
            println("| $(val) | — | — | — | — | — | $err |")
        else
            note = val >= 0.75 ? "hedge↓" : (val == 0.0 ? "max hedge" : "")
            xB   = isnan(row.s2.mean_xB) ? "?" : @sprintf("%.3f", row.s2.mean_xB)
            xBf  = isnan(row.s2.xB_frac) ? "?" : @sprintf("%.1f%%", row.s2.xB_frac * 100)
            @printf("| %.2f | %.2f | %.2f | %+.3f%% | %s | %s | %s |\n",
                    val, row.s1.V_rep, row.s2.V_rep, row.cev * 100, xB, xBf, note)
        end
    end
    println()
    println("_At rho_AB→1: cross-location returns perfectly correlated → no diversification value → hedge channel should collapse._")
end

function print_sweep_prelocate(dir)
    println("# p_relocate Sensitivity — v4 Option 1")
    println("Generated: $(Dates.now())")
    println()
    println("| p_reloc | V(E1_2L) | V(E2_2L) | CEV (%) | mean_xB@A | xB_frac | note |")
    println("|---------|----------|----------|---------|-----------|---------|------|")
    for (val, tag) in [(0.00,"0p00"), (0.06,"0p06"), (0.12,"0p12"), (0.30,"0p30")]
        row, err = cev_row_v4(dir, "E1_2L_preloc$(tag).json", "E2_2L_preloc$(tag).json")
        if row === nothing
            println("| $(val) | — | — | — | — | — | $err |")
        else
            note = val == 0.0 ? "no reloc→CEV??" : (val >= 0.20 ? "high mobility" : "")
            xB   = isnan(row.s2.mean_xB) ? "?" : @sprintf("%.3f", row.s2.mean_xB)
            xBf  = isnan(row.s2.xB_frac) ? "?" : @sprintf("%.1f%%", row.s2.xB_frac * 100)
            @printf("| %.2f | %.2f | %.2f | %+.3f%% | %s | %s | %s |\n",
                    val, row.s1.V_rep, row.s2.V_rep, row.cev * 100, xB, xBf, note)
        end
    end
    println()
    println("_At p_reloc=0: no relocation risk → hedge has no value; CEV should ≈ continuous-x channel only._")
end

function print_sweep_mortgage(dir)
    println("# Mortgage (LTV) Sensitivity — v4 Option 1")
    println("Generated: $(Dates.now())")
    println()
    println("| ltv_max | V(E1_2L) | V(E2_2L) | CEV (%) | mean_xB@A | note |")
    println("|---------|----------|----------|---------|-----------|------|")
    for (val, tag) in [(0.0,"0p0"), (0.5,"0p5"), (0.8,"0p8")]
        row, err = cev_row_v4(dir, "E1_2L_ltv$(tag).json", "E2_2L_ltv$(tag).json")
        if row === nothing
            println("| $(val) | — | — | — | — | $err |")
        else
            note = val >= 0.5 ? "YZ/Cocco range" : "no mortgage"
            xB   = isnan(row.s2.mean_xB) ? "?" : @sprintf("%.3f", row.s2.mean_xB)
            @printf("| %.1f | %.2f | %.2f | %+.3f%% | %s | %s |\n",
                    val, row.s1.V_rep, row.s2.V_rep, row.cev * 100, xB, note)
        end
    end
    println()
    println("_Higher LTV: mortgage substitutes for fractional x_ell → watch CEV shrink as LTV↑._")
end

function main()
    length(ARGS) < 2 && error("""
Usage:
  julia scripts/compute_cev_v4.jl baseline <e1.json> <e2.json>
  julia scripts/compute_cev_v4.jl rhoAB    <sweep_dir/>
  julia scripts/compute_cev_v4.jl prelocate <sweep_dir/>
  julia scripts/compute_cev_v4.jl mortgage  <sweep_dir/>
""")
    mode = ARGS[1]

    if mode == "baseline"
        length(ARGS) < 3 && error("baseline mode: provide <e1.json> <e2.json>")
        print_baseline_verdict(ARGS[2], ARGS[3])
    elseif mode == "rhoAB"
        print_sweep_rhoAB(ARGS[2])
    elseif mode == "prelocate"
        print_sweep_prelocate(ARGS[2])
    elseif mode == "mortgage"
        print_sweep_mortgage(ARGS[2])
    else
        error("Unknown mode '$mode'. Use: baseline, rhoAB, prelocate, mortgage")
    end
end

main()
