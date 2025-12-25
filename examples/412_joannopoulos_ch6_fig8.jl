# 3D Inverse Opal Photonic Crystal - Joannopoulos Chapter 6, Figure 8
# Air spheres in silicon (FCC lattice)
# Reference: Joannopoulos et al., "Photonic Crystals" 2nd ed., p.104
#
# Expected result: ~5% band gap between bands 8-9

using PhoXonic
using LinearAlgebra
using Plots
default(guidefontsize=14, tickfontsize=12, titlefontsize=14, left_margin=10Plots.mm, right_margin=10Plots.mm, top_margin=5Plots.mm, bottom_margin=10Plots.mm)

println("=" ^ 60)
println("Joannopoulos Ch.6 Fig.8: Inverse Opal")
println("=" ^ 60)
println()

# FCC lattice
a = 1.0
lat = fcc_lattice(a)

println("Structure: FCC lattice of air spheres in silicon")
println("  Primitive vectors:")
for (i, v) in enumerate(lat.vectors)
    println("    a$i = ", round.(v; digits=4))
end

# Joannopoulos parameters:
# - Background: silicon ε = 11.9
# - Air spheres: ε = 1
# - Radius: r ≈ 0.35a (close-packed)
# Note: Real inverse opals have interconnecting veins (omitted here)
silicon = Dielectric(11.9)
air = Dielectric(1.0)
radius = 0.35

geo = Geometry(lat, silicon, [(Sphere([0.0, 0.0, 0.0], radius), air)])

println("\nJoannopoulos parameters:")
println("  Background: silicon (ε = 11.9)")
println("  Air spheres: ε = 1, r/a = $radius")

# Solver setup - use TransverseEM (recommended)
resolution = (16, 16, 16)  # Match MPB resolution
cutoff = 7

println("\nSolver configuration:")
println("  Resolution: $resolution")
println("  Cutoff: $cutoff")
println("  Wave type: TransverseEM (2N basis)")

solver = Solver(TransverseEM(), geo, resolution, DenseMethod(); cutoff=cutoff)
println("  Plane waves: ", solver.basis.num_pw)
println("  Matrix size: ", 2 * solver.basis.num_pw, " × ", 2 * solver.basis.num_pw)

# K-path (Joannopoulos style)
kpath = kpath_fcc_joannopoulos(; a=a, npoints=15)
println("\nK-path: X → U|K → Γ → X → W → K")
println("K-points: ", length(kpath))

# Compute more bands (gap is between bands 8-9)
nbands = 12
println("\nComputing $nbands bands...")
println("(Gap expected between bands 8-9)")
bands_result = compute_bands(solver, kpath; bands=1:nbands, verbose=true)

# Find band gaps
println("\n=== Band Gaps ===")
gaps = find_all_gaps(bands_result; threshold=0.001)
if isempty(gaps)
    println("No complete band gap found.")
    println("Note: Inverse opal gap (~5%) requires high resolution")
else
    for g in gaps
        gap_pct = round(g.gap_ratio * 100; digits=1)
        println("Gap between bands $(g.bands): $(gap_pct)% gap-to-midgap")
    end
end

println("\n=== Expected (Joannopoulos Fig.8) ===")
println("~5% gap-to-midgap between bands 8-9")
println("Note: Our model omits interconnecting veins")

# Normalize frequencies
freqs_normalized = bands_result.frequencies ./ (2π)
dists = bands_result.distances

# Find band gap between bands 8-9 (if exists)
gap_bottom = maximum(freqs_normalized[:, 8])
gap_top = minimum(freqs_normalized[:, 9])
has_gap = gap_top > gap_bottom

# Manual plot with normalized frequencies (Joannopoulos style)
p = plot(;
    size=(700, 500),
    xlabel="Wave vector",
    ylabel="Frequency (ωa/2πc)",
    title="Inverse Opal - Joannopoulos Ch.6 Fig.8",
    legend=false,
    ylims=(0, 1.0),
    grid=true,
    framestyle=:box,
)

# Add band gap highlight (blue rectangle) if gap exists
if has_gap
    plot!(
        p,
        [dists[1], dists[end], dists[end], dists[1], dists[1]],
        [gap_bottom, gap_bottom, gap_top, gap_top, gap_bottom];
        fillrange=gap_bottom,
        fillalpha=0.3,
        fillcolor=:blue,
        linecolor=:transparent,
        linewidth=0,
        label="",
    )
end

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

output_file = joinpath(@__DIR__, "412_joannopoulos_ch6_fig8.png")
savefig(p, output_file)
println("\nSaved: $output_file")
