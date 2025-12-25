# MPB Diamond Lattice - Dielectric spheres in air
# Reference: https://mpb.readthedocs.io/en/latest/Data_Analysis_Tutorial/#diamond-lattice-of-spheres
#
# Expected result: 10.6% band gap between bands 2-3

using PhoXonic
using LinearAlgebra
using Plots
default(guidefontsize=14, tickfontsize=12, titlefontsize=14, left_margin=10Plots.mm, right_margin=10Plots.mm, top_margin=5Plots.mm, bottom_margin=10Plots.mm)

println("=" ^ 60)
println("MPB Diamond Lattice: Dielectric Spheres in Air")
println("=" ^ 60)
println()

# Diamond lattice = FCC with 2-atom basis
a = 1.0
lat = fcc_lattice(a)

println("Structure: Diamond lattice (FCC + 2-atom basis)")
println("  Primitive vectors:")
for (i, v) in enumerate(lat.vectors)
    println("    a$i = ", round.(v; digits=4))
end

# MPB parameters:
# - Background: air ε = 1
# - Dielectric spheres: ε = 11.56
# - Radius: r = 0.25a
air = Dielectric(1.0)
dielectric = Dielectric(11.56)
radius = 0.25

# Diamond basis: 2 spheres per primitive cell
# Positions in Cartesian coordinates (conventional cell units)
pos1 = [0.0, 0.0, 0.0]
pos2 = [a / 4, a / 4, a / 4]

geo = Geometry(
    lat, air, [(Sphere(pos1, radius), dielectric), (Sphere(pos2, radius), dielectric)]
)

println("\nMPB parameters:")
println("  Background: air (ε = 1)")
println("  Dielectric spheres: ε = 11.56, r/a = $radius")
println("  Basis: (0,0,0) and (a/4,a/4,a/4)")

# Solver setup - TransverseEM (no spurious longitudinal modes!)
resolution = (16, 16, 16)  # Match MPB resolution
cutoff = 7

println("\nSolver configuration:")
println("  Resolution: $resolution")
println("  Cutoff: $cutoff")
println("  Wave type: TransverseEM (2N basis, no spurious modes)")

solver = Solver(TransverseEM(), geo, resolution, DenseMethod(); cutoff=cutoff)
println("  Plane waves: ", solver.basis.num_pw)
println("  Matrix size: ", 2 * solver.basis.num_pw, " × ", 2 * solver.basis.num_pw)

# K-path: Γ → X → W → L → Γ → K
kpath = simple_kpath_fcc(; a=a, npoints=15)
println("\nK-path: Γ → X → W → L → Γ → K")
println("K-points: ", length(kpath))

# Compute bands - TransverseEM gives physical bands directly (no longitudinal modes)
nbands = 6
println("\nComputing bands 1:$nbands (TransverseEM: no spurious modes to skip!)...")
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

println("\n=== Expected (MPB tutorial) ===")
println("10.6% gap-to-midgap between bands 2-3")

# Manual plot (plot_bands has issues with band indexing)
freqs_norm = bands_result.frequencies ./ (2π)
dists = bands_result.distances

p = plot(;
    size=(700, 500),
    xlabel="Wave vector",
    ylabel="Frequency ωa/2πc",
    title="MPB Diamond: Dielectric Spheres (ε=11.56, r=0.25a)",
    legend=false,
)

# Highlight band gaps
for g in gaps
    gmin = g.max_lower / (2π)
    gmax = g.min_upper / (2π)
    hspan!(p, [gmin, gmax]; alpha=0.2, color=:yellow, label="")
end

# Plot each band
for b in 1:size(freqs_norm, 2)
    scatter!(p, dists, freqs_norm[:, b]; markersize=3, color=:blue)
end

# Add vertical lines at high-symmetry points
for (idx, label) in bands_result.labels
    vline!(p, [dists[idx]]; color=:gray, linestyle=:dash, alpha=0.5)
end

# Add x-axis labels
label_positions = [dists[idx] for (idx, _) in bands_result.labels]
label_names = [label for (_, label) in bands_result.labels]
plot!(p; xticks=(label_positions, label_names))

output_file = joinpath(@__DIR__, "413_mpb_diamond.png")
savefig(p, output_file)
println("\nSaved: $output_file")
println(
    "Frequency range: $(round(minimum(freqs_norm); digits=3)) - $(round(maximum(freqs_norm); digits=3))",
)
