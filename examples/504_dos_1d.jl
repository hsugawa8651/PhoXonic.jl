# 1D Photonic Crystal: Band Structure and DOS
# Shows band gaps and corresponding DOS features
#
# Note: RSKGF is currently implemented for 2D/3D.
# This example demonstrates basic DOS calculation for 1D.

using PhoXonic
using Plots

default(
    guidefontsize=14,
    tickfontsize=12,
    titlefontsize=14,
    left_margin=10Plots.mm,
    right_margin=10Plots.mm,
)

println("Example 504: 1D Photonic Crystal DOS")
println("=" ^ 50)

# 1D photonic crystal - Bragg stack
lat = lattice_1d(1.0)
geo = Geometry(lat, Dielectric(1.0), [
    (Segment(0.0, 0.3), Dielectric(9.0))
])
solver = Solver(Photonic1D(), geo, (64,); cutoff=10)

println("Structure: 1D Bragg stack")
println("Solver: Photonic1D, 64 grid, cutoff=10")
println("Plane waves: $(solver.basis.num_pw)")

# Frequency points for DOS
ω_values = range(0.1, 5.0, length=100)

# K-points for DOS sampling
k_points = collect(range(0.0, 0.5, length=20))

println("Frequencies: $(length(ω_values)) points")
println("K-points: $(length(k_points))")
println()

# Compute DOS
println("Computing DOS...")
t_dos = @elapsed dos_result = compute_dos(solver, ω_values, k_points; η=0.05)
println("  Time: $(round(t_dos, digits=2)) s")

# Compute band structure (1D uses direct Vector of k-points)
println("Computing band structure...")
kpath = collect(range(0.0, 0.5, length=41))
nbands = 5
nk = length(kpath)
frequencies = zeros(Float64, nk, nbands)

for (ik, k) in enumerate(kpath)
    ω, _ = solve(solver, k; bands=1:nbands)
    frequencies[ik, :] = ω
end

# Find band gaps manually
gaps = []
for b in 1:(nbands-1)
    max_lower = maximum(frequencies[:, b])
    min_upper = minimum(frequencies[:, b+1])
    gap = max(0.0, min_upper - max_lower)
    if gap > 0.01
        midgap = (min_upper + max_lower) / 2
        gap_ratio = gap / midgap
        push!(gaps, (bands=(b, b+1), gap=gap, gap_ratio=gap_ratio, min_upper=min_upper, max_lower=max_lower))
    end
end

println("\nBand gaps found: $(length(gaps))")
for g in gaps
    gap_pct = round(g.gap_ratio * 100; digits=1)
    println("  Bands $(g.bands): $(gap_pct)% gap-to-midgap at ω = $(round((g.min_upper + g.max_lower)/2, digits=2))")
end

# Plot
p = plot(layout=(1, 2), size=(1000, 450))

# Left: Band structure
dists = kpath

for b in 1:size(frequencies, 2)
    plot!(p, dists, frequencies[:, b];
        subplot=1,
        label="",
        linewidth=1.5,
        color=:blue,
    )
end

# Shade band gap regions
for g in gaps
    gap_bottom = g.max_lower
    gap_top = g.min_upper
    plot!(p,
        [dists[1], dists[end], dists[end], dists[1], dists[1]],
        [gap_bottom, gap_bottom, gap_top, gap_top, gap_bottom];
        subplot=1,
        fillrange=gap_bottom,
        fillalpha=0.3,
        fillcolor=:lightgreen,
        linecolor=:green,
        linewidth=1,
        label="",
    )
end

# High-symmetry point labels for 1D
label_positions = [0.0, 0.5]
label_names = ["Γ", "X"]

# Add vertical lines at high-symmetry points
for pos in label_positions
    vline!(p, [pos]; subplot=1, color=:gray, linestyle=:dot, alpha=0.5, label="")
end

plot!(p;
    subplot=1,
    xlabel="Wave vector k (π/a)",
    ylabel="Frequency ω",
    title="1D Photonic Crystal Bands",
    xticks=(label_positions, label_names),
    legend=false,
)

# Right: DOS
plot!(p, dos_result, ω_values;
    subplot=2,
    label="",
    linewidth=2,
    color=:blue,
)

# Mark band gap regions in DOS plot
for g in gaps
    gap_bottom = g.max_lower
    gap_top = g.min_upper
    hspan!(p, [gap_bottom, gap_top];
        subplot=2,
        alpha=0.3,
        color=:lightgreen,
        label="",
    )
end

plot!(p;
    subplot=2,
    xlabel="DOS (arb. units)",
    ylabel="Frequency ω",
    title="Density of States",
    legend=false,
)

output_file = joinpath(@__DIR__, "504_dos_1d.png")
savefig(p, output_file)
println("\nSaved: $output_file")
