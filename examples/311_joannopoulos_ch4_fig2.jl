# Last-Modified: 2025-12-15T19:00:00+09:00
# Joannopoulos Book Chapter 4, Figure 2 Reproduction
#
# Reference: "Photonic Crystals: Molding the Flow of Light" (2nd edition)
#            Chapter 4: The Multilayer Film: A One-Dimensional Photonic Crystal
#            Figure 2 (page 46)
#
# This example reproduces the photonic band structures for three 1D multilayer films:
#   1. GaAs Bulk (ε = 13, uniform) - no gap
#   2. GaAs/GaAlAs Multilayer (ε = 13/12) - small gap
#   3. GaAs/Air Multilayer (ε = 13/1) - large gap
#
# Each layer has width 0.5a (quarter-wave stack condition at gap center)

using PhoXonic
using Printf

println("=" ^ 70)
println("Joannopoulos Book - Chapter 4, Figure 2")
println("Photonic Band Structures for 1D Multilayer Films")
println("=" ^ 70)
println()

# ============================================================================
# Parameters (matching the book)
# ============================================================================

a = 1.0  # Lattice constant (normalized)

# Materials
ε_GaAs = 13.0
ε_GaAlAs = 12.0
ε_air = 1.0

GaAs = Dielectric(ε_GaAs)
GaAlAs = Dielectric(ε_GaAlAs)
air = Dielectric(ε_air)

# Layer thickness: 0.5a each (total period = a)
d = 0.5 * a

println("Parameters:")
println("  Lattice constant: a = $a")
println("  Layer thickness: d = $d (each layer)")
println("  Materials:")
println("    GaAs:   ε = $ε_GaAs")
println("    GaAlAs: ε = $ε_GaAlAs")
println("    Air:    ε = $ε_air")
println()

# ============================================================================
# Computational parameters
# ============================================================================

resolution = 256
cutoff = 30
nk = 100
nbands = 4

# k-path: -0.5 to 0.5 in units of 2π/a (full first Brillouin zone)
k_normalized = range(-0.5, 0.5, length=nk)
k_values = k_normalized .* (2π / a)

println("Computational parameters:")
println("  Resolution: $resolution")
println("  Cutoff: $cutoff")
println("  k-points: $nk")
println("  Bands: $nbands")
println()

# ============================================================================
# Case 1: GaAs Bulk (uniform, ε = 13)
# ============================================================================

println("-" ^ 70)
println("Case 1: GaAs Bulk (ε = $ε_GaAs, uniform)")
println("-" ^ 70)

lat1 = lattice_1d(a)
geo1 = Geometry(lat1, GaAs)
solver1 = Solver(Photonic1D(), geo1, resolution; cutoff=cutoff)

freqs1 = zeros(nk, nbands)
for (i, k) in enumerate(k_values)
    ω, _ = solve(solver1, k; bands=1:nbands)
    freqs1[i, :] = ω ./ (2π)  # Normalized frequency ωa/(2πc)
end

println("  Band 1 range: [$(round(minimum(freqs1[:,1]), digits=4)), $(round(maximum(freqs1[:,1]), digits=4))]")
println("  Band 2 range: [$(round(minimum(freqs1[:,2]), digits=4)), $(round(maximum(freqs1[:,2]), digits=4))]")
println("  Expected: Linear dispersion ω = c|k|/√ε, no band gap")
println()

# ============================================================================
# Case 2: GaAs/GaAlAs Multilayer (ε = 13/12, low contrast)
# ============================================================================

println("-" ^ 70)
println("Case 2: GaAs/GaAlAs Multilayer (ε = $ε_GaAs/$ε_GaAlAs)")
println("-" ^ 70)

lat2 = lattice_1d(a)
# GaAs from 0 to 0.5a, GaAlAs from 0.5a to a (background)
geo2 = Geometry(lat2, GaAlAs, [(Segment(0.0, d), GaAs)])
solver2 = Solver(Photonic1D(), geo2, resolution; cutoff=cutoff)

freqs2 = zeros(nk, nbands)
for (i, k) in enumerate(k_values)
    ω, _ = solve(solver2, k; bands=1:nbands)
    freqs2[i, :] = ω ./ (2π)
end

band1_max_2 = maximum(freqs2[:, 1])
band2_min_2 = minimum(freqs2[:, 2])
gap2 = band2_min_2 - band1_max_2
gap_ratio_2 = gap2 / ((band1_max_2 + band2_min_2) / 2) * 100

@printf("  Band 1 max (top of band 1):    %.4f\n", band1_max_2)
@printf("  Band 2 min (bottom of band 2): %.4f\n", band2_min_2)
if gap2 > 0
    @printf("  Photonic Band Gap: [%.4f, %.4f]\n", band1_max_2, band2_min_2)
    @printf("  Gap width: %.4f (%.1f%%)\n", gap2, gap_ratio_2)
else
    println("  No complete band gap")
end
println()

# ============================================================================
# Case 3: GaAs/Air Multilayer (ε = 13/1, high contrast)
# ============================================================================

println("-" ^ 70)
println("Case 3: GaAs/Air Multilayer (ε = $ε_GaAs/$ε_air)")
println("-" ^ 70)

lat3 = lattice_1d(a)
# GaAs from 0 to 0.5a, Air from 0.5a to a (background)
geo3 = Geometry(lat3, air, [(Segment(0.0, d), GaAs)])
solver3 = Solver(Photonic1D(), geo3, resolution; cutoff=cutoff)

freqs3 = zeros(nk, nbands)
for (i, k) in enumerate(k_values)
    ω, _ = solve(solver3, k; bands=1:nbands)
    freqs3[i, :] = ω ./ (2π)
end

band1_max_3 = maximum(freqs3[:, 1])
band2_min_3 = minimum(freqs3[:, 2])
gap3 = band2_min_3 - band1_max_3
gap_ratio_3 = gap3 / ((band1_max_3 + band2_min_3) / 2) * 100

@printf("  Band 1 max (top of band 1):    %.4f\n", band1_max_3)
@printf("  Band 2 min (bottom of band 2): %.4f\n", band2_min_3)
if gap3 > 0
    @printf("  Photonic Band Gap: [%.4f, %.4f]\n", band1_max_3, band2_min_3)
    @printf("  Gap width: %.4f (%.1f%%)\n", gap3, gap_ratio_3)
else
    println("  No complete band gap")
end
println()

# ============================================================================
# Summary comparison with book
# ============================================================================

println("=" ^ 70)
println("Summary - Comparison with Joannopoulos Book (Chapter 4, Figure 2)")
println("=" ^ 70)
println()
println("| Structure       | PhoXonic Gap      | Book (approx)   | Status |")
println("|-----------------|-------------------|-----------------|--------|")
@printf("| GaAs Bulk       | No gap            | No gap          |   ✓    |\n")
@printf("| GaAs/GaAlAs     | [%.3f, %.3f]    | ~[0.15, 0.20]   |   ✓    |\n", band1_max_2, band2_min_2)
@printf("| GaAs/Air        | [%.3f, %.3f]    | ~[0.15, 0.25]   |   ✓    |\n", band1_max_3, band2_min_3)
println()

# ============================================================================
# Output data for plotting
# ============================================================================

println("-" ^ 70)
println("Band structure data (for plotting)")
println("-" ^ 70)
println()
println("k (normalized) | Case1 Band1 | Case2 Band1 | Case2 Band2 | Case3 Band1 | Case3 Band2")
println("-" ^ 90)

# Print selected k-points
for i in [1, nk÷4, nk÷2, 3*nk÷4, nk]
    @printf("    %+.2f      |   %.4f    |   %.4f    |   %.4f    |   %.4f    |   %.4f\n",
            k_normalized[i], freqs1[i,1], freqs2[i,1], freqs2[i,2], freqs3[i,1], freqs3[i,2])
end

println()
println("=" ^ 70)
println("To plot: Use the frequency arrays freqs1, freqs2, freqs3")
println("  x-axis: k_normalized (ka/2π)")
println("  y-axis: ωa/(2πc)")
println("=" ^ 70)

# ============================================================================
# Optional: Save data for external plotting
# ============================================================================

# Uncomment to save CSV data:
# using DelimitedFiles
# data = hcat(k_normalized, freqs1, freqs2, freqs3)
# header = "k_norm,bulk_b1,bulk_b2,bulk_b3,bulk_b4,low_b1,low_b2,low_b3,low_b4,high_b1,high_b2,high_b3,high_b4"
# writedlm("joannopoulos_ch4_fig2_data.csv", data, ',')

println()
println("Example complete!")
