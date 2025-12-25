# 2D Phononic Crystal - Kushwaha et al. (1993) Figure 1
# Ni cylinders in Al matrix, SH wave (transverse)
# Reference: Phys. Rev. Lett. 71, 2022 (1993)
#
# Expected result: Band gap between bands 1-2

using PhoXonic
using LinearAlgebra
using Plots
default(guidefontsize=14, tickfontsize=12, titlefontsize=14, left_margin=10Plots.mm, right_margin=10Plots.mm, top_margin=5Plots.mm, bottom_margin=10Plots.mm)

println("=" ^ 60)
println("Kushwaha et al. (1993) Figure 1: Ni/Al Phononic Crystal")
println("=" ^ 60)
println()

# Square lattice
a = 1.0  # normalized lattice constant
lat = square_lattice(a)

println("Structure: 2D square lattice")
println("  Lattice constant: a = $a")

# Material parameters from Kushwaha (1993)
# Ni: ρ = 8936 kg/m³, c44 = 7.54×10¹⁰ Pa
# Al: ρ = 2697 kg/m³, c44 = 2.79×10¹⁰ Pa
#
# For SH wave, only shear modulus μ = c44 matters
# We use normalized units: ρ₀ = ρ_Al, μ₀ = c44_Al

ρ_Ni = 8936.0 / 2697.0   # normalized by Al density
μ_Ni = 7.54e10 / 2.79e10  # normalized by Al shear modulus

ρ_Al = 1.0   # reference
μ_Al = 1.0   # reference

# IsotropicElastic needs λ and μ, but for SH wave only μ matters
# Set λ = 0 for simplicity (doesn't affect SH mode)
Ni = IsotropicElastic(; ρ=ρ_Ni, λ=0.0, μ=μ_Ni)
Al = IsotropicElastic(; ρ=ρ_Al, λ=0.0, μ=μ_Al)

println("\nMaterial parameters (normalized):")
println("  Ni cylinder: ρ = $(round(ρ_Ni, digits=3)), μ = $(round(μ_Ni, digits=3))")
println("  Al matrix:   ρ = $(round(ρ_Al, digits=3)), μ = $(round(μ_Al, digits=3))")

# Filling fraction f = 0.1
# f = π r² / a² → r = a √(f/π)
f = 0.1
radius = a * sqrt(f / π)

println("\nGeometry:")
println("  Filling fraction: f = $f")
println("  Cylinder radius: r/a = $(round(radius/a, digits=4))")

# Ni cylinders in Al matrix
geo = Geometry(lat, Al, [(Circle([0.0, 0.0], radius), Ni)])

# Solver setup - SH wave (transverse mode)
resolution = (32, 32)
cutoff = 10

println("\nSolver configuration:")
println("  Wave type: SHWave (transverse, out-of-plane)")
println("  Resolution: $resolution")
println("  Cutoff: $cutoff")

solver = Solver(SHWave(), geo, resolution, DenseMethod(); cutoff=cutoff)
println("  Plane waves: ", solver.basis.num_pw)

# K-path: Γ → X → M → Γ
kpath = simple_kpath_square(; a=a, npoints=30)
println("\nK-path: Γ → X → M → Γ")
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

println("\n=== Expected (Kushwaha Fig.1) ===")
println("Band gap between bands 1-2")

# Normalize frequencies: ω → ωa/2πc_t where c_t = √(μ_Al/ρ_Al) = 1
# Since our units give c_t = 1, normalized freq = ω/(2π)
freqs_normalized = bands_result.frequencies ./ (2π)
dists = bands_result.distances

# Find band gap for highlighting (if exists)
gap_bottom = maximum(freqs_normalized[:, 1])
gap_top = minimum(freqs_normalized[:, 2])
has_gap = gap_top > gap_bottom

# Plot (Kushwaha style)
p = plot(;
    size=(700, 500),
    xlabel="Wave vector",
    ylabel="Frequency (ωa/2πcₜ)",
    title="Kushwaha (1993) Fig.1: Ni cylinders in Al, f=0.1",
    legend=false,
    ylims=(0, 1.5),
    grid=true,
    framestyle=:box,
)

# Add band gap highlight if exists
if has_gap
    plot!(
        p,
        [dists[1], dists[end], dists[end], dists[1], dists[1]],
        [gap_bottom, gap_bottom, gap_top, gap_top, gap_bottom];
        fillrange=gap_bottom,
        fillalpha=0.3,
        fillcolor=:lightblue,
        linecolor=:transparent,
        linewidth=0,
        label="",
    )
end

# Plot bands as lines
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

output_file = joinpath(@__DIR__, "211_kushwaha1993_fig1.png")
savefig(p, output_file)
println("\nSaved: $output_file")
