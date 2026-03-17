# DOS calculation for 2D photonic crystal
# Compares band structure with density of states
#
# Note: For RSKGF/MatrixFreeGF methods, load ReducedShiftedKrylov.jl

using PhoXonic
using Plots

default(
    guidefontsize=14,
    tickfontsize=12,
    titlefontsize=14,
    left_margin=10Plots.mm,
    right_margin=10Plots.mm,
    bottom_margin=10Plots.mm,
)

println("Example 502: DOS Calculation")
println("=" ^ 50)

# 2D photonic crystal - square lattice with dielectric rods
lat = square_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [
    (Circle([0.5, 0.5], 0.3), Dielectric(9.0))
])
solver = Solver(TMWave(), geo, (32, 32); cutoff=5)

println("Structure: Square lattice with dielectric rods")
println("Solver: TMWave, 32x32 grid, cutoff=5")
println("Plane waves: $(solver.basis.num_pw)")

# Compute band structure
println("\nComputing band structure...")
kpath = simple_kpath_square(; npoints=31)
bands_result = compute_bands(solver, kpath; bands=1:6)

# Frequency range for DOS
ω_max = maximum(bands_result.frequencies) * 1.1
ω_values = range(0.05, ω_max, length=100)

# K-points for BZ sampling
k_points = [[i/8, j/8] for i in 0:7 for j in 0:7]

println("Computing DOS ($(length(ω_values)) frequencies, $(length(k_points)) k-points)...")
t_dos = @elapsed dos_result = compute_dos(solver, ω_values, k_points, DirectGF(); η=0.02)
println("  Time: $(round(t_dos, digits=2)) s")

# Find band gaps
gaps = find_all_gaps(bands_result; threshold=0.01)
println("\nBand gaps: $(length(gaps))")
for g in gaps
    gap_pct = round(g.gap_ratio * 100; digits=1)
    println("  Bands $(g.bands): $(gap_pct)% gap-to-midgap")
end

# Plot
p = plot(layout=(1, 2), size=(1000, 450))

# Left: Band structure
dists = bands_result.distances
freqs = bands_result.frequencies

for b in 1:size(freqs, 2)
    plot!(p, dists, freqs[:, b];
        subplot=1,
        label="",
        linewidth=1.5,
        color=:blue,
    )
end

# Shade band gaps
for g in gaps
    plot!(p,
        [dists[1], dists[end], dists[end], dists[1], dists[1]],
        [g.max_lower, g.max_lower, g.min_upper, g.min_upper, g.max_lower];
        subplot=1,
        fillrange=g.max_lower,
        fillalpha=0.3,
        fillcolor=:lightgreen,
        linecolor=:green,
        linewidth=1,
        label="",
    )
end

# High-symmetry point labels
label_positions = [dists[idx] for (idx, _) in bands_result.labels]
label_names = [label for (_, label) in bands_result.labels]
for pos in label_positions
    vline!(p, [pos]; subplot=1, color=:gray, linestyle=:dot, alpha=0.5, label="")
end

plot!(p;
    subplot=1,
    xlabel="Wave vector",
    ylabel="Frequency ω",
    title="Band Structure",
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

# Mark band gaps in DOS
for g in gaps
    hspan!(p, [g.max_lower, g.min_upper];
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

output_file = joinpath(@__DIR__, "502_dos_rsk.png")
savefig(p, output_file)
println("\nSaved: $output_file")
