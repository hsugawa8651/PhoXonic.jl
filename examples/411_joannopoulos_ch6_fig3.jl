# 3D Diamond Photonic Crystal - Joannopoulos Chapter 6, Figure 3
# Air spheres in dielectric, diamond lattice (FCC + 2-atom basis)
# Reference: Joannopoulos et al., "Photonic Crystals" 2nd ed., p.98
#
# Expected result: ~29% band gap between bands 2-3

using PhoXonic
using LinearAlgebra
using Plots
default(guidefontsize=14, tickfontsize=12, titlefontsize=14, left_margin=10Plots.mm, right_margin=10Plots.mm, top_margin=5Plots.mm, bottom_margin=10Plots.mm)

println("=" ^ 60)
println("Joannopoulos Ch.6 Fig.3: Diamond Lattice Air Spheres")
println("=" ^ 60)
println()

# Diamond lattice = FCC with 2-atom basis
# Conventional cell: a = 1
# Primitive vectors (FCC): a1=(0,a/2,a/2), a2=(a/2,0,a/2), a3=(a/2,a/2,0)
# Basis: atoms at (0,0,0) and (a/4,a/4,a/4)
a = 1.0
lat = fcc_lattice(a)

println("Structure: Diamond lattice")
println("  FCC primitive vectors:")
for (i, v) in enumerate(lat.vectors)
    println("    a$i = ", round.(v; digits=4))
end

# Joannopoulos parameters:
# - Background: dielectric ε = 13
# - Air spheres: ε = 1
# - Radius: r = 0.325a (optimized for max gap)
# - Spheres overlap to form connected air network
dielectric = Dielectric(13.0)
air = Dielectric(1.0)
radius = 0.325

# Two air spheres per primitive cell (diamond basis)
pos1 = [0.0, 0.0, 0.0]
pos2 = [a / 4, a / 4, a / 4]

geo = Geometry(lat, dielectric, [(Sphere(pos1, radius), air), (Sphere(pos2, radius), air)])

println("\nJoannopoulos parameters:")
println("  Background: ε = 13")
println("  Air spheres: ε = 1, r/a = $radius")
println("  Basis: (0,0,0) and (a/4,a/4,a/4)")

# Solver setup - use TransverseEM (recommended for 3D photonic crystals)
resolution = (16, 16, 16)  # Match MPB resolution
cutoff = 7

println("\nSolver configuration:")
println("  Resolution: $resolution")
println("  Cutoff: $cutoff")
println("  Wave type: TransverseEM (2N basis, no spurious modes)")

solver = Solver(TransverseEM(), geo, resolution, DenseMethod(); cutoff=cutoff)
println("  Plane waves: ", solver.basis.num_pw)
println("  Matrix size: ", 2 * solver.basis.num_pw, " × ", 2 * solver.basis.num_pw)

# K-path: Joannopoulos uses X → U|K → Γ → X → W → K
# Note: U and K are equivalent points in the FCC BZ
kpath = kpath_fcc_joannopoulos(; a=a, npoints=15)
println("\nK-path: X → U|K → Γ → X → W → K")
println("K-points: ", length(kpath))

# Compute bands
nbands = 6
println("\nComputing $nbands bands...")
bands_result = compute_bands(solver, kpath; bands=1:nbands, verbose=true)

# Find band gaps
println("\n=== Band Gaps ===")
gaps = find_all_gaps(bands_result; threshold=0.001)
if isempty(gaps)
    println("No complete band gap found.")
else
    for g in gaps
        gap_pct = round(g.gap_ratio * 100; digits=1)
        println("Gap between bands $(g.bands): $(gap_pct)% gap-to-midgap")
    end
end

println("\n=== Expected (Joannopoulos Fig.3) ===")
println("~29% gap-to-midgap between bands 2-3")

# Normalize frequencies: PhoXonic outputs ω, Joannopoulos uses ωa/2πc
# Since a=1, we need ω/(2π)
freqs_normalized = bands_result.frequencies ./ (2π)
dists = bands_result.distances

# Find band gap for highlighting
gap_bottom = maximum(freqs_normalized[:, 2])
gap_top = minimum(freqs_normalized[:, 3])

# Manual plot with normalized frequencies (Joannopoulos style)
p = plot(;
    size=(700, 500),
    xlabel="Wave vector",
    ylabel="Frequency (ωa/2πc)",
    title="Diamond Air Spheres - Joannopoulos Ch.6 Fig.3",
    legend=false,
    ylims=(0, 0.8),
    grid=true,
    framestyle=:box,
)

# Add band gap highlight (yellow rectangle)
plot!(
    p,
    [dists[1], dists[end], dists[end], dists[1], dists[1]],
    [gap_bottom, gap_bottom, gap_top, gap_top, gap_bottom];
    fillrange=gap_bottom,
    fillalpha=0.3,
    fillcolor=:yellow,
    linecolor=:transparent,
    linewidth=0,
    label="",
)

# Plot each band as lines (Joannopoulos style)
for b in 1:size(freqs_normalized, 2)
    plot!(p, dists, freqs_normalized[:, b]; linewidth=1.5, color=:blue, label="")
end

# Add vertical lines at high-symmetry points
for (idx, label) in bands_result.labels
    vline!(p, [dists[idx]]; color=:gray, linestyle=:solid, alpha=0.5, linewidth=0.5)
end

# Add x-axis labels
label_positions = [dists[idx] for (idx, _) in bands_result.labels]
label_names = [label for (_, label) in bands_result.labels]
plot!(p; xticks=(label_positions, label_names))

output_file = joinpath(@__DIR__, "411_joannopoulos_ch6_fig3.png")
savefig(p, output_file)
println("\nSaved: $output_file")
