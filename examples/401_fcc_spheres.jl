# 3D FCC Photonic Crystal Band Structure
# Dielectric spheres in air - TransverseEM (recommended for 3D)

using PhoXonic
using LinearAlgebra
using Plots
default(guidefontsize=14, tickfontsize=12, titlefontsize=14, left_margin=10Plots.mm, right_margin=10Plots.mm, top_margin=5Plots.mm, bottom_margin=10Plots.mm)

println("=== 3D FCC Photonic Crystal ===")
println("Dielectric spheres (ε=12) in air, r=0.25a")
println()

# FCC lattice
a = 1.0
lat = fcc_lattice(a)
println("FCC lattice (conventional cell a = $a):")
println("  Primitive vectors:")
for (i, v) in enumerate(lat.vectors)
    println("    a$i = ", round.(v; digits=4))
end

# Geometry: air background with dielectric sphere at origin
air = Dielectric(1.0)
sphere_material = Dielectric(12.0)
radius = 0.25
geo = Geometry(lat, air, [(Sphere([0.0, 0.0, 0.0], radius), sphere_material)])

println("\nGeometry: sphere with r/a = $radius")

# Create solver with TransverseEM (recommended for 3D photonic crystals)
# TransverseEM uses 2N×2N matrices with no spurious longitudinal modes
println("\nCreating TransverseEM solver...")
resolution = (16, 16, 16)  # Match MPB resolution
cutoff = 5

solver = Solver(TransverseEM(), geo, resolution, DenseMethod(); cutoff=cutoff)
println("  Resolution: $resolution")
println("  Cutoff: $cutoff")
println("  Plane waves: ", solver.basis.num_pw)
println("  Matrix size: ", 2 * solver.basis.num_pw, " × ", 2 * solver.basis.num_pw)

# Create k-path: Γ → X → W → L → Γ → K
println("\nComputing band structure...")
kpath = simple_kpath_fcc(; a=a, npoints=15)
println("K-path: Γ → X → W → L → Γ → K")
println("Number of k-points: ", length(kpath))

# Compute bands
nbands = 6
println("\nComputing $nbands bands...")
bands_result = compute_bands(solver, kpath; bands=1:nbands, verbose=true)
println("Computation complete.")

# Find band gaps
println("\n=== Band Gaps ===")
gaps = find_all_gaps(bands_result; threshold=0.001)
if isempty(gaps)
    println("No significant band gaps found.")
else
    for g in gaps
        println(
            "Gap between bands ",
            g.bands,
            ": Δω = ",
            round(g.gap; digits=4),
            " (",
            round(g.gap_ratio * 100; digits=1),
            "% gap-to-midgap)",
        )
    end
end

# Print label positions
println("\n=== High Symmetry Points ===")
dists = bands_result.distances
for (idx, label) in bands_result.labels
    println(label, " at distance ", round(dists[idx]; digits=4))
end

# Plot band structure
println("\n=== Generating Plot ===")
p = plot_bands(
    bands_result;
    title="3D FCC Photonic Crystal (ε=12 spheres, r=0.25a)",
    ylabel="Frequency ω (2πc/a)",
    scatter=true,
    markersize=2,
    color=:blue,
    size=(700, 500),
)

# Save figure
output_file = joinpath(@__DIR__, "401_fcc_bands.png")
savefig(p, output_file)
println("Saved: $output_file")

println("\n✅ 3D FCC photonic crystal example complete!")
