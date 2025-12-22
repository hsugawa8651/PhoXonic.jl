# Last-Modified: 2025-12-15T19:00:00+09:00
# Defect Mode in a 2D Photonic Crystal
# Demonstrates LDOS calculation to find localized defect modes
#
# A point defect (missing rod) in a square lattice of dielectric rods
# creates localized modes in the photonic bandgap.

using PhoXonic
using LinearAlgebra
using StaticArrays
using Plots

println("=== Defect Mode in 2D Photonic Crystal ===")
println("Supercell approach with LDOS calculation")
println()

# ============================================================================
# Perfect Crystal (for reference bandgap)
# ============================================================================

println("Step 1: Computing perfect crystal band structure...")

a = 1.0  # Lattice constant
lat = square_lattice(a)

# Materials: Silicon rods in air
air = Dielectric(1.0)
silicon = Dielectric(11.7)

# Rod radius (filling fraction ~ 0.28)
r = 0.3

# Perfect crystal geometry
geo_perfect = Geometry(lat, air, [(Circle([0.0, 0.0], r), silicon)])

# Solver for perfect crystal (TM mode has larger gaps for rods in air)
solver_perfect = Solver(TMWave(), geo_perfect, (32, 32); cutoff=5)

# Compute band structure
kpath = simple_kpath_square(; a=a, npoints=30)
bands_perfect = compute_bands(solver_perfect, kpath; bands=1:8)

# Find bandgaps
gaps = find_all_gaps(bands_perfect; threshold=0.01)
println("\nPerfect crystal bandgaps (TM):")
for g in gaps
    println(
        "  Bands $(g.bands): gap = $(round(g.gap, digits=4)), " *
        "range = [$(round(g.max_lower, digits=4)), $(round(g.min_upper, digits=4))]",
    )
end

# Select the first significant gap for defect mode analysis
if !isempty(gaps)
    gap = gaps[1]
    ω_min = gap.max_lower
    ω_max = gap.min_upper
    ω_mid = (ω_min + ω_max) / 2
    println("\nUsing gap between bands $(gap.bands)")
    println("  Gap range: [$ω_min, $ω_max]")
    println("  Gap midpoint: $ω_mid")
else
    # Fallback if no gap found
    ω_mid = 0.35 * 2π / a
    ω_min = 0.3 * 2π / a
    ω_max = 0.4 * 2π / a
    println("\nNo significant gap found, using fallback range")
end

# ============================================================================
# Supercell with Point Defect
# ============================================================================

println("\nStep 2: Creating supercell with point defect...")

# Supercell size (NxN unit cells)
N_super = 5

# Create supercell lattice
lat_super = Lattice(SVector(N_super * a, 0.0), SVector(0.0, N_super * a))

# Create inclusions: rods at all positions except center
inclusions = Tuple{Circle,Dielectric}[]
center_idx = (N_super + 1) / 2

for i in 1:N_super
    for j in 1:N_super
        # Skip center rod (point defect)
        if i == ceil(Int, center_idx) && j == ceil(Int, center_idx)
            continue
        end
        # Position of this rod
        x = (i - 0.5) * a
        y = (j - 0.5) * a
        push!(inclusions, (Circle([x, y], r), silicon))
    end
end

geo_defect = Geometry(lat_super, air, inclusions)

println("  Supercell: $(N_super) x $(N_super) unit cells")
println("  Missing rod at center (defect)")
println("  Total rods: $(length(inclusions)) / $(N_super^2)")

# ============================================================================
# LDOS at Defect Site
# ============================================================================

println("\nStep 3: Computing LDOS at defect site...")

# Solver for supercell
# Use lower resolution per unit cell to keep computation manageable
resolution_per_cell = 16
resolution = (N_super * resolution_per_cell, N_super * resolution_per_cell)
solver_defect = Solver(TMWave(), geo_defect, resolution; cutoff=5)

println("  Resolution: $resolution")
println("  Plane waves: $(solver_defect.basis.num_pw)")

# Defect position (center of supercell)
defect_pos = [N_super * a / 2, N_super * a / 2]

# Frequency range for LDOS (focusing on the gap)
n_freq = 50
ω_range = range(ω_min * 0.8, ω_max * 1.2; length=n_freq)

# K-points for Brillouin zone sampling (reduced zone due to supercell)
# For supercell, Γ point is often sufficient
k_points = [[0.0, 0.0]]

# Broadening parameter
η = 0.01

println("  Frequency range: [$(minimum(ω_range)), $(maximum(ω_range))]")
println("  Computing LDOS...")

# Use the standard (dense) LDOS for this example
# For larger systems, use compute_ldos(..., MatrixFreeGF())
ldos = compute_ldos(solver_defect, defect_pos, collect(ω_range), k_points; η=η)

# Find peak (defect mode)
ldos_max, idx_max = findmax(ldos)
ω_defect = ω_range[idx_max]

println("\nResults:")
println("  Maximum LDOS at ω = $(round(ω_defect, digits=4))")
println("  This is the defect mode frequency")

# Check if it's in the gap
if ω_min < ω_defect < ω_max
    println("  Defect mode is INSIDE the bandgap!")
else
    println("  Note: Defect mode may be outside the gap region")
end

# ============================================================================
# Visualization
# ============================================================================

println("\nStep 4: Creating plots...")

# Plot 1: Perfect crystal band structure
p_bands = plot(;
    xlabel="Wave vector",
    ylabel="Frequency (2πc/a)",
    title="Perfect Crystal TM Bands",
    legend=false,
    grid=true,
    size=(500, 400),
)

dists = bands_perfect.distances
for b in 1:size(bands_perfect.frequencies, 2)
    plot!(
        p_bands, dists, bands_perfect.frequencies[:, b] * a / (2π); linewidth=2, color=:blue
    )
end

# Add gap region
if !isempty(gaps)
    hspan!(
        p_bands,
        [ω_min * a / (2π), ω_max * a / (2π)];
        alpha=0.2,
        color=:yellow,
        label="Bandgap",
    )
end

# Add labels
label_positions = [dists[i] for (i, _) in bands_perfect.labels]
label_names = [l for (_, l) in bands_perfect.labels]
vline!(p_bands, label_positions; color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_bands, label_positions, label_names)

# Plot 2: LDOS
p_ldos = plot(
    collect(ω_range) * a / (2π),
    ldos;
    xlabel="Frequency (2πc/a)",
    ylabel="LDOS (arb. units)",
    title="LDOS at Defect Site",
    linewidth=2,
    color=:red,
    legend=false,
    grid=true,
    size=(500, 400),
)

# Mark the gap region
if !isempty(gaps)
    vspan!(p_ldos, [ω_min * a / (2π), ω_max * a / (2π)]; alpha=0.2, color=:yellow)
end

# Mark the defect mode
vline!(p_ldos, [ω_defect * a / (2π)]; color=:green, linestyle=:dash, linewidth=2)
annotate!(
    p_ldos, ω_defect * a / (2π), ldos_max * 0.9, text("Defect\nmode", :left, 8, :green)
)

# Combined plot
p_combined = plot(p_bands, p_ldos; layout=(1, 2), size=(1000, 400))

# Save plots
savefig(p_bands, joinpath(@__DIR__, "501_defect_bands.png"))
savefig(p_ldos, joinpath(@__DIR__, "501_defect_ldos.png"))
savefig(p_combined, joinpath(@__DIR__, "501_defect_combined.png"))

println("\nSaved: 501_defect_bands.png")
println("Saved: 501_defect_ldos.png")
println("Saved: 501_defect_combined.png")

# ============================================================================
# Comparison: LDOS inside vs outside defect
# ============================================================================

println("\nStep 5: Comparing LDOS at different positions...")

# Position away from defect (in a rod)
rod_pos = [0.5 * a, 0.5 * a]

ldos_rod = compute_ldos(solver_defect, rod_pos, collect(ω_range), k_points; η=η)

p_comparison = plot(;
    xlabel="Frequency (2πc/a)",
    ylabel="LDOS (arb. units)",
    title="LDOS Comparison: Defect vs Bulk",
    legend=:topright,
    grid=true,
    size=(600, 400),
)

plot!(
    p_comparison,
    collect(ω_range) * a / (2π),
    ldos;
    label="At defect",
    linewidth=2,
    color=:red,
)
plot!(
    p_comparison,
    collect(ω_range) * a / (2π),
    ldos_rod;
    label="In bulk",
    linewidth=2,
    color=:blue,
)

if !isempty(gaps)
    vspan!(
        p_comparison,
        [ω_min * a / (2π), ω_max * a / (2π)];
        alpha=0.2,
        color=:yellow,
        label="Bandgap",
    )
end

savefig(p_comparison, joinpath(@__DIR__, "501_defect_ldos_comparison.png"))
println("Saved: 501_defect_ldos_comparison.png")

display(p_combined)

println("\n=== Done ===")
