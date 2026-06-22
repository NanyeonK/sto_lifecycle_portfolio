#!/usr/bin/env julia
# Export lifecycle policy means and V-slice from serialized v4 solver results.
#
# Produces (per regime):
#   output/diagnostics/p6_option1_e{1,2}_policy_summary.csv
#   output/diagnostics/p6_option1_e{1,2}_v_slice_t1.csv
#
# Run from repo root after server1 Gate 1 baselines:
#   julia scripts/export_policy_csv.jl
#
# SAVE_PATH must have been set when running run_option1_e{1,2}.sh
# so that p6_option1_e{1,2}_result.jls exist.

using Serialization, Statistics

OUTDIR = "output/diagnostics"
LOC_A  = 1
LOC_B  = 2

function export_policy_csv(result_path::String, policy_out::String, vslice_out::String)
    if !isfile(result_path)
        println("SKIP (missing): $result_path — run server1 baseline with SAVE_PATH set")
        return
    end
    println("Loading: $result_path")
    result = open(deserialize, result_path)

    T    = size(result.value, 1) - 1   # skip terminal pseudo-period
    ages = collect(25:(25 + T - 1))

    # ── Policy summary: mean x_A, x_B by age and location ──────────────────
    rows = Tuple[]
    for t in 1:T
        for (iell, loc_lbl) in [(LOC_A, "ellA"), (LOC_B, "ellB")]
            # Entry state (x_prev = 0 for both dims, index 1)
            f_e   = view(result.feasible,  t, :, :, iell, 1, 1)
            xAp_e = view(result.xA_policy, t, :, :, iell, 1, 1)
            xBp_e = view(result.xB_policy, t, :, :, iell, 1, 1)
            n_e   = count(f_e)

            # All x_prev states averaged
            f_all  = result.feasible[ t, :, :, iell, :, :]
            xAp_a  = result.xA_policy[t, :, :, iell, :, :][f_all]
            xBp_a  = result.xB_policy[t, :, :, iell, :, :][f_all]
            n_all  = count(f_all)

            push!(rows, (
                ages[t], loc_lbl,
                n_e  > 0 ? mean(xAp_e[f_e]) : NaN,
                n_e  > 0 ? mean(xBp_e[f_e]) : NaN,
                n_all > 0 ? mean(xAp_a) : NaN,
                n_all > 0 ? mean(xBp_a) : NaN,
                n_all > 0 ? count(xBp_a .> 0.0) / n_all : NaN,
                n_all,
            ))
        end
    end

    open(policy_out, "w") do io
        println(io, "age,loc,mean_xA_entry,mean_xB_entry,mean_xA_all,mean_xB_all,frac_xB_pos,n_feasible")
        for r in rows
            println(io, join(r, ","))
        end
    end
    println("  Policy summary → $policy_out  ($(length(rows)) rows)")

    # ── V-slice at t=1, all w, iz_mid, ell=A, x_prev=(0,0) ─────────────────
    nw    = size(result.value, 2)
    iz_mid = max(1, div(size(result.value, 3), 2))

    vslice_e0 = fill(NaN, nw)
    vslice_e1 = fill(NaN, nw)
    vslice_e2 = fill(NaN, nw)
    feasible_flag = fill(false, nw)

    # Use E2_2L regime value if available; ell=A, iz_mid, x_prev=(0,0)
    for iw in 1:nw
        feasible_flag[iw] = result.feasible[1, iw, iz_mid, LOC_A, 1, 1]
        vslice_e2[iw] = result.value[1, iw, iz_mid, LOC_A, 1, 1]
    end

    open(vslice_out, "w") do io
        println(io, "iw,v_t1_ellA,feasible")
        for iw in 1:nw
            println(io, "$iw,$(vslice_e2[iw]),$(feasible_flag[iw])")
        end
    end
    println("  V-slice (t=1, ellA, iz_mid, xprev=0) → $vslice_out  ($nw rows)")
end

export_policy_csv(
    "$OUTDIR/p6_option1_e2_result.jls",
    "$OUTDIR/p6_option1_e2_policy_summary.csv",
    "$OUTDIR/p6_option1_e2_v_slice_t1.csv",
)
export_policy_csv(
    "$OUTDIR/p6_option1_e1_result.jls",
    "$OUTDIR/p6_option1_e1_policy_summary.csv",
    "$OUTDIR/p6_option1_e1_v_slice_t1.csv",
)

println("\nExport complete. Next:")
println("  python scripts/fig5_mean_x_age.py")
println("  python scripts/plot_lifecycle_profiles.py")
println("  python scripts/fig4_v_slice.py")
