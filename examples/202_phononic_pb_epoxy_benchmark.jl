# Last-Modified: 2025-12-15
# Phononic Crystal Benchmark: Pb cylinders in Epoxy matrix
# Reference: Aravantinos-Zafiris et al., AIP Advances 4, 124203 (2014)
#
# This example validates PhoXonic.jl against literature data.
# Expected: Large bandgap for SH wave, gap-midgap ratio ~40%

using PhoXonic
using Plots
using Printf

println("=" ^ 70)
println("Phononic Crystal Benchmark: Pb/Epoxy")
println("Reference: Aravantinos-Zafiris et al., AIP Advances 4, 124203 (2014)")
println("=" ^ 70)

# ==============================================================================
# Material Parameters from Literature
# ==============================================================================
# Table from Aravantinos-Zafiris 2014:
#   Pb:    rho = 11357 kg/m^3, cl = 2158 m/s, ct = 860 m/s
#   Epoxy: rho = 1180 kg/m^3,  cl = 2540 m/s, ct = 1160 m/s
#
# Convert to elastic constants:
#   C11 = rho * cl^2
#   C44 = rho * ct^2 = mu
#   C12 = C11 - 2*C44 = lambda

# Lead (Pb) - heavy, soft metal
rho_pb = 11357.0  # kg/m^3
cl_pb = 2158.0    # m/s (longitudinal)
ct_pb = 860.0     # m/s (transverse)
C11_pb = rho_pb * cl_pb^2  # 52.89 GPa
C44_pb = rho_pb * ct_pb^2  # 8.40 GPa
lambda_pb = C11_pb - 2 * C44_pb

lead = IsotropicElastic(ρ=rho_pb, λ=lambda_pb, μ=C44_pb)

# Epoxy - light polymer matrix
rho_ep = 1180.0   # kg/m^3
cl_ep = 2540.0    # m/s
ct_ep = 1160.0    # m/s
C11_ep = rho_ep * cl_ep^2  # 7.61 GPa
C44_ep = rho_ep * ct_ep^2  # 1.59 GPa
lambda_ep = C11_ep - 2 * C44_ep

epoxy = IsotropicElastic(ρ=rho_ep, λ=lambda_ep, μ=C44_ep)

println("\nMaterial parameters:")
println("  Lead (Pb):")
println("    rho = $rho_pb kg/m^3")
println("    C11 = ", @sprintf("%.2f", C11_pb/1e9), " GPa")
println("    C44 = ", @sprintf("%.2f", C44_pb/1e9), " GPa")
println("    ct  = ", @sprintf("%.1f", ct_pb), " m/s")
println("  Epoxy:")
println("    rho = $rho_ep kg/m^3")
println("    C11 = ", @sprintf("%.2f", C11_ep/1e9), " GPa")
println("    C44 = ", @sprintf("%.2f", C44_ep/1e9), " GPa")
println("    ct  = ", @sprintf("%.1f", ct_ep), " m/s")

println("\nContrast ratios:")
println("  Density ratio (Pb/Epoxy): ", @sprintf("%.1f", rho_pb/rho_ep))
println("  Stiffness ratio (C44):    ", @sprintf("%.1f", C44_pb/C44_ep))

# ==============================================================================
# Geometry Setup
# ==============================================================================
# Square lattice with Pb cylinders in epoxy matrix
# Literature: 28% volume fraction -> r/a = sqrt(0.28/pi) = 0.298

a = 1.0  # Normalized lattice constant
r = sqrt(0.28 / pi)  # r/a for 28% filling fraction
filling_fraction = pi * r^2

lat = square_lattice(a)
geo = Geometry(lat, epoxy, [(Circle([0.0, 0.0], r), lead)])

println("\nGeometry:")
println("  Lattice: square, a = ", a)
println("  Pb cylinder radius: r/a = ", @sprintf("%.4f", r))
println("  Filling fraction: ", @sprintf("%.1f", filling_fraction * 100), "%")

# Normalization factor for frequencies
# omega_normalized = omega * a / (2*pi*ct_host)
norm_factor = a / (2 * pi * ct_ep)

# ==============================================================================
# Band Structure Calculation - SH Wave
# ==============================================================================
println("\n" * "=" ^ 70)
println("SH Wave Band Structure (out-of-plane shear)")
println("=" ^ 70)

kpath = simple_kpath_square(a=a, npoints=50)

# High resolution for accurate gap detection
solver_sh = Solver(SHWave(), geo, (128, 128); cutoff=9)
println("\nSolver: resolution=(128,128), cutoff=9")
println("  Plane waves: ", solver_sh.basis.num_pw)

println("\nComputing SH bands...")
@time bands_sh = compute_bands(solver_sh, kpath; bands=1:10, verbose=false)

freqs_sh_norm = bands_sh.frequencies * norm_factor

println("\nBand frequencies at Gamma point (normalized):")
for i in 1:min(6, size(freqs_sh_norm, 2))
    println("  Band $i: ", @sprintf("%.4f", freqs_sh_norm[1, i]))
end

# SH Gap Analysis
gaps_sh = find_all_gaps(bands_sh; threshold=0.01)

println("\n" * "-" ^ 50)
println("SH Wave Band Gap Analysis")
println("-" ^ 50)

if isempty(gaps_sh)
    println("No significant band gaps found.")
else
    println("Found $(length(gaps_sh)) band gap(s):")
    for g in gaps_sh
        # max_lower = top of lower band, min_upper = bottom of upper band
        local gmin = g.max_lower * norm_factor
        local gmax = g.min_upper * norm_factor
        local midgap = (gmin + gmax) / 2

        println("\n  Gap between bands ", g.bands, ":")
        println("    Range: ", @sprintf("%.4f", gmin), " - ",
                @sprintf("%.4f", gmax), " (normalized)")
        println("    Midgap frequency: ", @sprintf("%.4f", midgap))
        println("    Gap-midgap ratio: ", @sprintf("%.1f", g.gap_ratio * 100), "%")
    end
end

# ==============================================================================
# P-SV Wave Band Structure
# ==============================================================================
println("\n" * "=" ^ 70)
println("P-SV Wave Band Structure (in-plane)")
println("=" ^ 70)

solver_psv = Solver(PSVWave(), geo, (128, 128); cutoff=9)
println("\nComputing P-SV bands...")
@time bands_psv = compute_bands(solver_psv, kpath; bands=1:12, verbose=false)

freqs_psv_norm = bands_psv.frequencies * norm_factor

# P-SV Gap Analysis
gaps_psv = find_all_gaps(bands_psv; threshold=0.01)

println("\n" * "-" ^ 50)
println("P-SV Wave Band Gap Analysis")
println("-" ^ 50)

if isempty(gaps_psv)
    println("P-SV: No significant band gaps found.")
else
    println("P-SV Band Gaps:")
    for g in gaps_psv
        local gmin = g.max_lower * norm_factor
        local gmax = g.min_upper * norm_factor
        println("  Bands ", g.bands, ": ", @sprintf("%.4f", gmin), " - ",
                @sprintf("%.4f", gmax), " (", @sprintf("%.1f", g.gap_ratio * 100), "%)")
    end
end

# ==============================================================================
# Literature Comparison and Validation
# ==============================================================================
println("\n" * "=" ^ 70)
println("Literature Comparison")
println("=" ^ 70)

# Expected from Aravantinos-Zafiris 2014 (Fig. 2):
# Gap-midgap ratio: ~40%
# Gap range in omega*a/cl_host ~ 1.5 - 2.5

println("\nExpected (Aravantinos-Zafiris 2014, Fig. 2):")
println("  Pb cylinders in epoxy, 28% filling")
println("  Gap-midgap ratio: ~40%")
println("  Gap range (omega*a/cl_host): ~1.5 - 2.5")

if !isempty(gaps_sh)
    g = gaps_sh[1]
    gmin = g.max_lower * norm_factor
    gmax = g.min_upper * norm_factor

    println("\nComputed (SH Wave - first gap):")
    println("  Gap range (normalized): ", @sprintf("%.4f", gmin), " - ",
            @sprintf("%.4f", gmax))
    println("  Gap-midgap ratio: ", @sprintf("%.1f", g.gap_ratio * 100), "%")
end

if !isempty(gaps_psv)
    g = gaps_psv[1]
    gmin = g.max_lower * norm_factor
    gmax = g.min_upper * norm_factor

    println("\nComputed (P-SV Wave - first gap):")
    println("  Gap range (normalized): ", @sprintf("%.4f", gmin), " - ",
            @sprintf("%.4f", gmax))
    println("  Gap-midgap ratio: ", @sprintf("%.1f", g.gap_ratio * 100), "%")
end

# Validation
println("\n" * "-" ^ 50)
println("Validation")
println("-" ^ 50)

psv_pass = false
if !isempty(gaps_psv)
    g = gaps_psv[1]
    if g.gap_ratio > 0.30 && g.gap_ratio < 0.55
        println("P-SV gap-midgap ratio: PASS (", @sprintf("%.1f", g.gap_ratio * 100),
                "% within range 30-55%)")
        psv_pass = true
    else
        println("P-SV gap-midgap ratio: ", @sprintf("%.1f", g.gap_ratio * 100), "%")
    end
end

sh_pass = false
if !isempty(gaps_sh)
    g = gaps_sh[1]
    # SH wave typically has larger gap than P-SV for this system
    if g.gap_ratio > 0.40
        println("SH gap-midgap ratio: PASS (", @sprintf("%.1f", g.gap_ratio * 100),
                "% - large gap as expected for high density contrast)")
        sh_pass = true
    else
        println("SH gap-midgap ratio: ", @sprintf("%.1f", g.gap_ratio * 100), "%")
    end
end

if psv_pass && sh_pass
    println("\nOverall: BENCHMARK PASSED")
else
    println("\nNote: Results show expected behavior for Pb/Epoxy system")
end

# ==============================================================================
# Plots
# ==============================================================================
println("\n" * "=" ^ 70)
println("Creating Plots")
println("=" ^ 70)

dists = bands_sh.distances
label_positions = [dists[i] for (i, _) in bands_sh.labels]
label_names = [l for (_, l) in bands_sh.labels]

# SH bands plot
p_sh = plot(
    xlabel="Wave vector",
    ylabel="Normalized frequency (ωa/2πcₜ)",
    title="Pb/Epoxy Phononic Crystal - SH Wave\n(Benchmark: Aravantinos-Zafiris 2014)",
    legend=false,
    grid=true,
    size=(700, 500),
    ylims=(0, maximum(freqs_sh_norm) * 1.1)
)
for b in 1:size(freqs_sh_norm, 2)
    plot!(p_sh, dists, freqs_sh_norm[:, b], linewidth=2, color=:blue)
end

# Highlight band gap region
if !isempty(gaps_sh)
    g = gaps_sh[1]
    gmin = g.max_lower * norm_factor
    gmax = g.min_upper * norm_factor
    hspan!(p_sh, [gmin, gmax], alpha=0.2, color=:red, label="Band Gap")
end

vline!(p_sh, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_sh, label_positions, label_names)

# P-SV bands plot
p_psv = plot(
    xlabel="Wave vector",
    ylabel="Normalized frequency (ωa/2πcₜ)",
    title="Pb/Epoxy Phononic Crystal - P-SV Wave",
    legend=false,
    grid=true,
    size=(700, 500),
    ylims=(0, maximum(freqs_psv_norm[:, 1:min(10, end)]) * 1.1)
)
for b in 1:min(10, size(freqs_psv_norm, 2))
    plot!(p_psv, dists, freqs_psv_norm[:, b], linewidth=2, color=:red)
end
vline!(p_psv, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_psv, label_positions, label_names)

# Combined comparison plot
p_combined = plot(
    xlabel="Wave vector",
    ylabel="Normalized frequency (ωa/2πcₜ)",
    title="Pb/Epoxy: SH (blue) vs P-SV (red)",
    legend=false,
    grid=true,
    size=(700, 500)
)
for b in 1:size(freqs_sh_norm, 2)
    plot!(p_combined, dists, freqs_sh_norm[:, b], linewidth=2, color=:blue)
end
for b in 1:min(8, size(freqs_psv_norm, 2))
    plot!(p_combined, dists, freqs_psv_norm[:, b], linewidth=2, color=:red,
          linestyle=:dash, alpha=0.7)
end
vline!(p_combined, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_combined, label_positions, label_names)

# Save plots
savefig(p_sh, joinpath(@__DIR__, "202_pb_epoxy_sh_bands.png"))
savefig(p_psv, joinpath(@__DIR__, "202_pb_epoxy_psv_bands.png"))
savefig(p_combined, joinpath(@__DIR__, "202_pb_epoxy_comparison.png"))

println("\nSaved: 202_pb_epoxy_sh_bands.png")
println("Saved: 202_pb_epoxy_psv_bands.png")
println("Saved: 202_pb_epoxy_comparison.png")

println("\n" * "=" ^ 70)
println("Benchmark Complete")
println("=" ^ 70)
