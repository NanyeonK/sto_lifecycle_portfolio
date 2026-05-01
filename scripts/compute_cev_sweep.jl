#!/usr/bin/env julia
# Compute CEV sensitivity table from VFI summary JSON pairs.
#
# Usage:
#   julia scripts/compute_cev_sweep.jl <outdir> <sweep_type>
#
# sweep_type in {rhoAB, prelocate, txcost}
#
# Reads pairs of E1_2L and E2_2L summary JSONs from <outdir>,
# computes CEV at the midpoint state (iw_mid, iz_mid, ell=A),
# and writes a Markdown table to stdout.
#
# CEV formula: (V_E2 / V_E1)^(1/(1-γ)) − 1
# For γ=5:      (V_E2 / V_E1)^(−1/4) − 1
# Both V values are negative (CRRA with γ>1); V_E2 > V_E1 (less negative) implies ratio < 1,
# giving CEV > 0 after the negative exponent flips the inequality.

using JSON3, Printf

function compute_cev(V_E2::Float64, V_E1::Float64, gamma::Float64)::Float64
    ratio = V_E2 / V_E1
    ratio <= 0.0 && error("V_E2/V_E1 = $ratio <= 0 — check signs of value function outputs")
    return ratio^(1.0 / (1.0 - gamma)) - 1.0
end

function read_summary(path::String)
    d = JSON3.read(read(path, String))
    V  = Float64(d["V_t1_midpoint_ellA"])
    gm = Float64(d["params"]["gamma"])
    ts = Float64(d["params"]["tau_sell"])
    tb = Float64(d["params"]["tau_buy"])
    ap = Bool(d["params"]["apply_tau_buy_at_reloc"])
    rhoAB = Float64(d["params"]["rho_AB"])
    preloc = Float64(d["params"]["p_relocate_working"])
    return (; V, gm, ts, tb, ap, rhoAB, preloc)
end

function cev_row(dir, tag_e1, tag_e2)
    f_e1 = joinpath(dir, tag_e1)
    f_e2 = joinpath(dir, tag_e2)
    isfile(f_e1) || return nothing, "MISSING: $f_e1"
    isfile(f_e2) || return nothing, "MISSING: $f_e2"
    s1 = read_summary(f_e1)
    s2 = read_summary(f_e2)
    cev = compute_cev(s2.V, s1.V, s1.gm)
    return (; V_E1=s1.V, V_E2=s2.V, cev, gamma=s1.gm, s1.ts, s1.tb, s1.ap, s1.rhoAB, s1.preloc), nothing
end

function main()
    length(ARGS) < 2 && error("Usage: compute_cev_sweep.jl <outdir> <sweep_type>")
    outdir = ARGS[1]
    sweep  = ARGS[2]

    println("# CEV sensitivity: $sweep")
    println("Generated: $(Dates.now())")
    println()

    if sweep == "rhoAB"
        println("| rho_AB | V(E1_2L) | V(E2_2L) | CEV (%) | note |")
        println("|--------|----------|----------|---------|------|")
        for (val, tag) in [(0.00,"0p00"), (0.25,"0p25"), (0.50,"0p50"), (0.75,"0p75"), (0.95,"0p95")]
            row, err = cev_row(outdir, "E1_2L_rhoAB$(tag).json", "E2_2L_rhoAB$(tag).json")
            if row === nothing
                println("| $(val) | — | — | — | $err |")
            else
                note = val >= 0.75 ? "hedge↓" : (val == 0.0 ? "max hedge" : "")
                @printf("| %.2f | %.2f | %.2f | %+.3f%% | %s |\n",
                        val, row.V_E1, row.V_E2, row.cev * 100, note)
            end
        end

    elseif sweep == "prelocate"
        println("| p_reloc | V(E1_2L) | V(E2_2L) | CEV (%) | note |")
        println("|---------|----------|----------|---------|------|")
        for (val, tag) in [(0.00,"0p00"), (0.02,"0p02"), (0.06,"0p06"), (0.12,"0p12")]
            row, err = cev_row(outdir, "E1_2L_preloc$(tag).json", "E2_2L_preloc$(tag).json")
            if row === nothing
                println("| $(val) | — | — | — | $err |")
            else
                note = val == 0.0 ? "no reloc→CEV→0" : (val >= 0.10 ? "high mobility" : "")
                @printf("| %.2f | %.2f | %.2f | %+.3f%% | %s |\n",
                        val, row.V_E1, row.V_E2, row.cev * 100, note)
            end
        end

    elseif sweep == "txcost"
        println("| scenario | tau_sell | tau_buy | apply_buy | V(E1_2L) | V(E2_2L) | CEV (%) |")
        println("|----------|----------|---------|-----------|----------|----------|---------|")
        scenarios = [
            ("notx",   "NOTX (0+0)"),
            ("sell6",  "sell 6% only"),
            ("rt8p5",  "sell 6%+buy 2.5% = 8.5% rt"),
            ("rt10",   "sell 6%+buy 4% = 10% rt"),
            ("rt12",   "sell 6%+buy 6% = 12% rt"),
        ]
        for (tag, label) in scenarios
            row, err = cev_row(outdir, "E1_2L_$(tag).json", "E2_2L_$(tag).json")
            if row === nothing
                println("| $label | — | — | — | — | — | $err |")
            else
                @printf("| %s | %.3f | %.3f | %s | %.2f | %.2f | %+.3f%% |\n",
                        label, row.ts, row.tb, row.ap ? "yes" : "no",
                        row.V_E1, row.V_E2, row.cev * 100)
            end
        end
    else
        error("Unknown sweep type '$sweep'. Use: rhoAB, prelocate, txcost")
    end

    println()
    println("_CEV = ((V_E2/V_E1)^(1/(1−γ)) − 1)×100.  State: midpoint (iw_mid, iz_mid, ell=A)._")
    println("_tau_buy approximation: owner who relocates pays tau_buy deducted from relocation wealth._")
    println("_E2_2L: tokens portable across moves — no tau_sell, no tau_buy at relocation._")
end

using Dates
main()
