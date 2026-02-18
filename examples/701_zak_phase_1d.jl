# Example 701: 1D Zak Phase Calculation
#
# This example demonstrates calculating the Zak phase for a 1D photonic crystal.
# The Zak phase is quantized to 0 or π for systems with inversion symmetry.

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

# Create a 1D photonic crystal with a dielectric layer
# Unit cell: air-dielectric-air
lat = lattice_1d(1.0)

# Geometry A: dielectric layer centered in unit cell
geo_A = Geometry(lat, Dielectric(1.0), [(Segment(0.35, 0.65), Dielectric(9.0))])

# Geometry B: dielectric layer at edges (same material, different origin)
geo_B = Geometry(lat, Dielectric(1.0), [
    (Segment(0.0, 0.15), Dielectric(9.0)),
    (Segment(0.85, 1.0), Dielectric(9.0)),
])

# Create solvers
solver_A = Solver(Photonic1D(), geo_A, 128; cutoff = 20)
solver_B = Solver(Photonic1D(), geo_B, 128; cutoff = 20)

# Compute Zak phases for bands 1-4
println("Computing Zak phases...")
println()

n_bands = 4
zak_A = compute_zak_phase(solver_A, 1:n_bands; n_k = 100)
zak_B = compute_zak_phase(solver_B, 1:n_bands; n_k = 100)

println("Geometry A (centered dielectric):")
for (i, phase) in enumerate(zak_A.phases)
    phase_norm = phase / π
    println("  Band $i: Zak phase = $(round(phase, digits=4)) ($(round(phase_norm, digits=2))π)")
end
println()

println("Geometry B (edge dielectric):")
for (i, phase) in enumerate(zak_B.phases)
    phase_norm = phase / π
    println("  Band $i: Zak phase = $(round(phase, digits=4)) ($(round(phase_norm, digits=2))π)")
end
println()

# Plot Zak phases comparison
p = plot(
    title = "Zak Phase Comparison",
    xlabel = "Band index",
    ylabel = "Zak phase / π",
    legend = :topright,
    ylims = (-1.2, 1.2),
    yticks = ([-1, -0.5, 0, 0.5, 1], ["-π", "-π/2", "0", "π/2", "π"]),
    xticks = 1:n_bands,
    size = (700, 500),
)

# Add horizontal lines at 0 and ±π (behind)
hline!([0, 1, -1], color = :gray, linestyle = :dash, label = "")

# Plot Geometry B - red crosses (×)
scatter!(
    1:n_bands,
    zak_B.phases ./ π,
    label = "Geometry B (edge)",
    markersize = 10,
    markershape = :xcross,
    color = :red,
    markerstrokewidth = 2,
)

# Plot Geometry A - blue plus (+)
scatter!(
    1:n_bands,
    zak_A.phases ./ π,
    label = "Geometry A (centered)",
    markersize = 10,
    markershape = :cross,
    color = :blue,
    markerstrokewidth = 2,
)

savefig(p, "701_zak_phase_1d.png")
println("Saved: 701_zak_phase_1d.png")

# Also compute and plot band structure for reference
println()
println("Computing band structure...")

k_values = range(0.0, 0.5, length = 51)
n_k = length(k_values)
n_plot_bands = 6
freqs = zeros(n_k, n_plot_bands)

for (i, k) in enumerate(k_values)
    ω = solve_at_k(solver_A, k, solver_A.method; bands = 1:n_plot_bands)
    freqs[i, :] = ω
end

p2 = plot(
    title = "1D Photonic Crystal Band Structure",
    xlabel = "Wave vector k (π/a)",
    ylabel = "Frequency (ωa/2πc)",
    legend = false,
    size = (700, 500),
)

for band in 1:n_plot_bands
    plot!(p2, k_values .* 2, freqs[:, band], color = :blue, linewidth = 2)
end

savefig(p2, "701_zak_phase_1d_bands.png")
println("Saved: 701_zak_phase_1d_bands.png")

println()
println("Note: For inversion-symmetric systems, Zak phase is quantized to 0 or π.")
println("      Different unit cell choices yield different Zak phases.")
