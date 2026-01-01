# Example 702: 2D Wilson Loop Spectrum
#
# This example demonstrates calculating the Wilson loop spectrum for a 2D photonic crystal.
# The Wilson spectrum reveals topological properties through the winding of phases.

using PhoXonic
using Plots
default(
    guidefontsize = 14,
    tickfontsize = 12,
    titlefontsize = 14,
    legendfontsize = 11,
    left_margin = 10Plots.mm,
    right_margin = 10Plots.mm,
    top_margin = 5Plots.mm,
    bottom_margin = 10Plots.mm,
)

# Create a 2D photonic crystal: square lattice with circular rod
lat = square_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Circle([0.0, 0.0], 0.3), Dielectric(9.0))])

# Create solver for TM polarization
solver = Solver(TMWave(), geo, (32, 32); cutoff = 5)

println("Computing Wilson loop spectrum for 2D photonic crystal...")
println("Structure: square lattice, circular rod (r=0.3a), ε=9.0")
println()

# Compute Wilson spectrum for bands 1-2
# Loop direction :b2 means we scan along b1 and compute loops in b2 direction
result = compute_wilson_spectrum(solver, 1:2; n_k_path = 41, n_k_loop = 50, loop_direction = :b2)

println("Wilson spectrum calculated:")
println("  k-path points: $(length(result.k_values))")
println("  Number of bands: $(length(result.bands))")
println("  Loop direction: $(result.loop_direction)")
println()

# Calculate winding numbers
println("Winding numbers:")
for band in 1:length(result.bands)
    w = winding_number(result, band)
    println("  Band $band: winding = $w")
end
println()

# Plot Wilson spectrum
p1 = plot(
    title = "Wilson Loop Spectrum (TM mode, bands 1-2)",
    xlabel = "k₁ (2π/a)",
    ylabel = "Wilson phase / π",
    legend = :topright,
    ylims = (-1.1, 1.1),
    yticks = ([-1, -0.5, 0, 0.5, 1], ["-π", "-π/2", "0", "π/2", "π"]),
    size = (700, 500),
)

# Plot band 2 - red crosses (×)
scatter!(
    p1,
    result.k_values,
    result.phases[:, 2] ./ π,
    label = "Band 2",
    markersize = 8,
    markershape = :xcross,
    color = :red,
    markerstrokewidth = 2,
)

# Plot band 1 - blue plus (+)
scatter!(
    p1,
    result.k_values,
    result.phases[:, 1] ./ π,
    label = "Band 1",
    markersize = 8,
    markershape = :cross,
    color = :blue,
    markerstrokewidth = 2,
)

hline!([0], color = :gray, linestyle = :dash, label = "")

savefig(p1, "702_wilson_loop_2d.png")
println("Saved: 702_wilson_loop_2d.png")

# Also plot band structure for reference
println()
println("Computing band structure...")

kpath = kpath_square()
bands = compute_bands(solver, kpath; bands = 1:6)

# Manual plot with clear axis labels
dists = bands.distances
freqs = bands.frequencies
label_positions = Float64[dists[i] for (i, _) in bands.labels]
label_names = String[l for (_, l) in bands.labels]

p2 = plot(;
    xlabel = "Wave vector",
    ylabel = "Frequency (ωa/2πc)",
    title = "2D Photonic Crystal (TM mode)",
    legend = false,
    grid = true,
    size = (700, 500),
)
for b in 1:size(freqs, 2)
    plot!(p2, dists, freqs[:, b]; linewidth = 2, color = :blue)
end
vline!(p2, label_positions; color = :gray, linestyle = :dash, alpha = 0.5)
xticks!(p2, label_positions, label_names)

savefig(p2, "702_wilson_loop_2d_bands.png")
println("Saved: 702_wilson_loop_2d_bands.png")

println()
println("Note: Non-zero winding indicates non-trivial topology.")
println("      This simple structure has trivial topology (winding = 0).")
