# Example 604: PWE vs TMM Comparison
#
# Compares band structure calculations from Plane Wave Expansion (PWE)
# and Transfer Matrix Method (TMM) for 1D photonic and phononic crystals.
# Requires: Plots.jl
#
# NOTE: 1D phononic PWE shows larger errors due to missing inverse rule.
#       See ../TODO-1d-phononic-inverse-rule.md for details.
#       When inverse rule is implemented, verify this example converges.

using PhoXonic
using Plots

# ============================================================================
# Photonic Crystal: Dielectric Stack
# ============================================================================

println("=== Photonic Crystal: PWE vs TMM ===\n")

# Material parameters
n_hi = 3.0   # High index (e.g., Si)
n_lo = 1.5   # Low index (e.g., SiO2)

mat_hi = Dielectric(n_hi^2)
mat_lo = Dielectric(n_lo^2)

# Period and fill fraction
a = 1.0           # Normalized period
ff = 0.3          # Fill fraction of high-index material
d_hi = ff * a
d_lo = (1 - ff) * a

println("Structure: high-n ($n_hi) / low-n ($n_lo) dielectric stack")
println("Period: a = $a")
println("Fill fraction: $ff")
println()

# ----------------------------------------
# TMM setup
# ----------------------------------------
unit_cell_tmm = [Layer(mat_hi, d_hi), Layer(mat_lo, d_lo)]
ml = Multilayer(unit_cell_tmm, mat_lo, mat_lo)
solver_tmm = TMMSolver(Photonic1D(), ml)

# ----------------------------------------
# PWE setup
# ----------------------------------------
lat = lattice_1d(a)
geo = Geometry(lat, mat_lo, [(Segment(0.0, d_hi), mat_hi)])
solver_pwe = Solver(Photonic1D(), geo, 256; cutoff=20)

# ----------------------------------------
# Compute bands
# ----------------------------------------
n_bands = 5
k_points = 21

bands_tmm = tmm_bandstructure(solver_tmm; k_points=k_points, bands=1:n_bands)

# PWE: compute at same k-points
k_values = range(0, π/a, length=k_points)
ω_pwe = zeros(k_points, n_bands)

for (ik, k) in enumerate(k_values)
    ω, _ = solve(solver_pwe, k; bands=1:n_bands)
    ω_pwe[ik, :] = ω
end

# ----------------------------------------
# Compare results
# ----------------------------------------
println("Comparison at zone boundary (k = π/a):\n")
println("Band    TMM ω      PWE ω      Diff (%)")
println("-" ^ 40)

for i in 1:n_bands
    ω_tmm_edge = bands_tmm.frequencies[end, i]
    ω_pwe_edge = ω_pwe[end, i]
    diff_pct = 100 * abs(ω_tmm_edge - ω_pwe_edge) / ω_pwe_edge
    println("  $i    $(round(ω_tmm_edge, digits=4))    $(round(ω_pwe_edge, digits=4))    $(round(diff_pct, digits=2))%")
end

# ----------------------------------------
# Bandgap comparison
# ----------------------------------------
println("\nBandgap between bands 1-2:")

# TMM
gap1_max_tmm = maximum(bands_tmm.frequencies[:, 1])
gap2_min_tmm = minimum(bands_tmm.frequencies[:, 2])
gap_tmm = gap2_min_tmm - gap1_max_tmm

# PWE
gap1_max_pwe = maximum(ω_pwe[:, 1])
gap2_min_pwe = minimum(ω_pwe[:, 2])
gap_pwe = gap2_min_pwe - gap1_max_pwe

println("  TMM: gap = $(round(gap_tmm, digits=4)) ($(round(gap1_max_tmm, digits=4)) to $(round(gap2_min_tmm, digits=4)))")
println("  PWE: gap = $(round(gap_pwe, digits=4)) ($(round(gap1_max_pwe, digits=4)) to $(round(gap2_min_pwe, digits=4)))")

# ============================================================================
# Phononic Crystal: Steel/Epoxy
# ============================================================================

println("\n" * "=" ^ 60)
println("=== Phononic Crystal: PWE vs TMM ===\n")

# Materials
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)

# Period
a_phon = 0.01  # 10 mm
ff_phon = 0.25
d_steel = ff_phon * a_phon
d_epoxy = (1 - ff_phon) * a_phon

println("Structure: Steel/Epoxy multilayer")
println("Period: a = $(a_phon*1e3) mm")
println("Fill fraction (steel): $ff_phon")
println()

# ----------------------------------------
# TMM setup
# ----------------------------------------
unit_cell_phon = [Layer(steel, d_steel), Layer(epoxy, d_epoxy)]
ml_phon = Multilayer(unit_cell_phon, epoxy, epoxy)
solver_tmm_phon = TMMSolver(Longitudinal1D(), ml_phon)

# ----------------------------------------
# PWE setup
# ----------------------------------------
lat_phon = lattice_1d(a_phon)
geo_phon = Geometry(lat_phon, epoxy, [(Segment(0.0, d_steel), steel)])
solver_pwe_phon = Solver(Longitudinal1D(), geo_phon, 256; cutoff=20)

# ----------------------------------------
# Compute bands
# ----------------------------------------
n_bands_phon = 4

bands_tmm_phon = tmm_bandstructure(solver_tmm_phon; k_points=k_points, bands=1:n_bands_phon)

# PWE: compute at same k-points
k_values_phon = range(0, π/a_phon, length=k_points)
ω_pwe_phon = zeros(k_points, n_bands_phon)

for (ik, k) in enumerate(k_values_phon)
    ω, _ = solve(solver_pwe_phon, k; bands=1:n_bands_phon)
    ω_pwe_phon[ik, :] = ω
end

# ----------------------------------------
# Compare results
# ----------------------------------------
println("Comparison at zone boundary (k = π/a):\n")
println("Band    TMM ω (rad/s)    PWE ω (rad/s)    Diff (%)")
println("-" ^ 55)

for i in 1:n_bands_phon
    ω_tmm_edge = bands_tmm_phon.frequencies[end, i]
    ω_pwe_edge = ω_pwe_phon[end, i]
    diff_pct = 100 * abs(ω_tmm_edge - ω_pwe_edge) / ω_pwe_edge
    println("  $i    $(round(ω_tmm_edge, sigdigits=5))    $(round(ω_pwe_edge, sigdigits=5))    $(round(diff_pct, digits=1))%")
end

# ----------------------------------------
# Bandgap comparison
# ----------------------------------------
println("\nBandgap between bands 1-2:")

# TMM
gap1_max_tmm_phon = maximum(bands_tmm_phon.frequencies[:, 1])
gap2_min_tmm_phon = minimum(bands_tmm_phon.frequencies[:, 2])
gap_tmm_phon = gap2_min_tmm_phon - gap1_max_tmm_phon

# PWE
gap1_max_pwe_phon = maximum(ω_pwe_phon[:, 1])
gap2_min_pwe_phon = minimum(ω_pwe_phon[:, 2])
gap_pwe_phon = gap2_min_pwe_phon - gap1_max_pwe_phon

println("  TMM: gap = $(round(gap_tmm_phon, sigdigits=4)) rad/s")
println("  PWE: gap = $(round(gap_pwe_phon, sigdigits=4)) rad/s")

gap_center_tmm = (gap1_max_tmm_phon + gap2_min_tmm_phon) / 2
gap_center_pwe = (gap1_max_pwe_phon + gap2_min_pwe_phon) / 2
rel_gap_tmm = gap_tmm_phon / gap_center_tmm
rel_gap_pwe = gap_pwe_phon / gap_center_pwe

println("  TMM: relative gap = $(round(100*rel_gap_tmm, digits=1))%")
println("  PWE: relative gap = $(round(100*rel_gap_pwe, digits=1))%")

# ============================================================================
# Summary
# ============================================================================

println("\n" * "=" ^ 60)
println("=== Summary ===\n")

println("Photonic: PWE and TMM produce consistent band structures.")
println("  - PWE converges to TMM with increasing cutoff")
println("  - Differences < 0.1% with cutoff=100")
println()
println("Phononic: PWE shows larger differences from TMM.")
println("  - Current PWE uses direct rule (C11), not inverse rule (C11_inv)")
println("  - High contrast materials (Steel/Epoxy ~37:1) cause Gibbs phenomenon")
println("  - TODO: Implement inverse rule for 1D phononic PWE")
println()
println("TMM advantages:")
println("  - Exact for 1D multilayers")
println("  - Direct R/T spectrum calculation")
println("  - Handles lossy materials (complex ε)")
println()
println("PWE advantages:")
println("  - Generalizes to 2D/3D")
println("  - Direct band structure")
println("  - Handles arbitrary geometries")

# ============================================================================
# Plot 1: Photonic Band Structure Comparison
# ============================================================================

println("\nGenerating plots...")

p1 = plot(
    xlabel="Wave vector ka/π",
    ylabel="Frequency ωa/2πc",
    title="Photonic Crystal: PWE vs TMM",
    legend=:topleft,
    grid=true,
    size=(800, 500)
)

# Normalize k-values to [0, 1] (ka/π)
k_norm = collect(k_values) .* a ./ π

# Plot TMM bands
for b in 1:n_bands
    plot!(p1, k_norm, bands_tmm.frequencies[:, b],
          linewidth=2, color=:blue, label=(b==1 ? "TMM" : ""))
end

# Plot PWE bands
for b in 1:n_bands
    plot!(p1, k_norm, ω_pwe[:, b],
          linewidth=2, linestyle=:dash, color=:red, label=(b==1 ? "PWE" : ""))
end

# Highlight first bandgap
if gap_tmm > 0
    hspan!(p1, [gap1_max_tmm, gap2_min_tmm], alpha=0.15, color=:yellow, label="Gap")
end

# Add zone labels
xticks!(p1, [0, 1], ["Γ", "X"])

savefig(p1, joinpath(@__DIR__, "604_photonic_comparison.png"))
println("Saved: 604_photonic_comparison.png")

# ============================================================================
# Plot 2: Phononic Band Structure Comparison
# ============================================================================

p2 = plot(
    xlabel="Wave vector ka/π",
    ylabel="Frequency ω (rad/s)",
    title="Phononic Crystal: PWE vs TMM",
    legend=:topleft,
    grid=true,
    size=(800, 500)
)

# Normalize k-values
k_norm_phon = collect(k_values_phon) .* a_phon ./ π

# Plot TMM bands
for b in 1:n_bands_phon
    plot!(p2, k_norm_phon, bands_tmm_phon.frequencies[:, b],
          linewidth=2, color=:blue, label=(b==1 ? "TMM" : ""))
end

# Plot PWE bands
for b in 1:n_bands_phon
    plot!(p2, k_norm_phon, ω_pwe_phon[:, b],
          linewidth=2, linestyle=:dash, color=:red, label=(b==1 ? "PWE" : ""))
end

# Highlight first bandgap (TMM)
if gap_tmm_phon > 0
    hspan!(p2, [gap1_max_tmm_phon, gap2_min_tmm_phon], alpha=0.15, color=:yellow, label="Gap (TMM)")
end

# Add zone labels
xticks!(p2, [0, 1], ["Γ", "X"])

savefig(p2, joinpath(@__DIR__, "604_phononic_comparison.png"))
println("Saved: 604_phononic_comparison.png")

# ============================================================================
# Plot 3: Error Analysis (Photonic)
# ============================================================================

p3 = plot(
    xlabel="Band number",
    ylabel="Relative difference (%)",
    title="PWE vs TMM: Relative Error at Zone Boundary",
    legend=:topleft,
    grid=true,
    size=(800, 500)
)

# Photonic errors
errors_photonic = Float64[]
for i in 1:n_bands
    ω_tmm_edge = bands_tmm.frequencies[end, i]
    ω_pwe_edge = ω_pwe[end, i]
    push!(errors_photonic, 100 * abs(ω_tmm_edge - ω_pwe_edge) / ω_pwe_edge)
end

# Phononic errors
errors_phononic = Float64[]
for i in 1:n_bands_phon
    ω_tmm_edge = bands_tmm_phon.frequencies[end, i]
    ω_pwe_edge = ω_pwe_phon[end, i]
    push!(errors_phononic, 100 * abs(ω_tmm_edge - ω_pwe_edge) / ω_pwe_edge)
end

bar!(p3, 1:n_bands, errors_photonic, label="Photonic", alpha=0.7, color=:blue)
bar!(p3, (1:n_bands_phon) .+ 0.3, errors_phononic, label="Phononic", alpha=0.7, color=:red, bar_width=0.3)

savefig(p3, joinpath(@__DIR__, "604_error_analysis.png"))
println("Saved: 604_error_analysis.png")

display(p1)
println("\nDone!")
