# Last-Modified: 2025-12-15T19:00:00+09:00
# 3D FCC Photonic Crystal Band Structure
# Dielectric spheres in air - Full vector EM (H-field formulation)

using PhoXonic
using LinearAlgebra
using Plots

println("=== 3D FCC Photonic Crystal ===")
println("Dielectric spheres (ε=12) in air, r=0.25a")
println()

# FCC lattice
a = 1.0
lat = fcc_lattice(a)
println("FCC lattice (conventional cell a = $a):")
println("  Primitive vectors:")
for (i, v) in enumerate(lat.vectors)
    println("    a$i = ", round.(v, digits=4))
end

# Geometry: air background with dielectric sphere at origin
air = Dielectric(1.0)
sphere_material = Dielectric(12.0)  # High dielectric like GaAs
radius = 0.25  # in units of a
geo = Geometry(lat, air, [(Sphere([0.0, 0.0, 0.0], radius), sphere_material)])

println("\nGeometry: sphere with r/a = $radius")

# Create solver with shift-and-invert (required for 3D to skip longitudinal modes)
println("\nCreating FullVectorEM solver...")
resolution = (12, 12, 12)
cutoff = 7
shift = 0.01

solver = Solver(FullVectorEM(), geo, resolution, KrylovKitMethod(shift=shift); cutoff=cutoff)
println("  Resolution: $resolution")
println("  Cutoff: $cutoff")
println("  Number of plane waves: ", solver.basis.num_pw)
println("  Total DOF: ", solver.basis.num_pw * 3, " (3 components)")
println("  Shift (σ): $shift")

# Create k-path: Γ → X → W → L → Γ → K (simplified connected path)
# Note: At Γ point (k=0), the lowest transverse modes also have ω→0,
# which can cause issues with shift-and-invert. We use a small offset.
println("\nComputing band structure...")
kpath = simple_kpath_fcc(a=a, npoints=20)
println("K-path: Γ → X → W → L → Γ → K")
println("Number of k-points: ", length(kpath))
println("Note: Γ point may show anomalous values due to k=0 singularity")

# Compute bands (fewer bands for 3D due to computational cost)
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
        println("Gap between bands ", g.bands, ": Δω = ", round(g.gap, digits=4),
                " (", round(g.gap_ratio*100, digits=1), "% gap-to-midgap)")
    end
end

# Output data for plotting
println("\n=== Band Structure Data ===")
println("# distance  ", join(["band$i" for i in 1:nbands], "  "))
freqs = bands_result.frequencies
dists = bands_result.distances
for i in 1:length(dists)
    print(round(dists[i], digits=4))
    for b in 1:nbands
        print("  ", round(freqs[i, b], digits=4))
    end
    println()
end

# Print label positions
println("\n=== High Symmetry Points ===")
for (idx, label) in bands_result.labels
    println(label, " at distance ", round(dists[idx], digits=4))
end

# Comparison with different solver methods (optional)
println("\n=== Solver Method Comparison ===")
k_test = [0.5, 0.5, 0.0]  # Test point
println("Testing at k = $k_test")

# KrylovKit (shift-and-invert finds eigenvalues closest to σ)
solver_krylov = Solver(FullVectorEM(), geo, resolution, KrylovKitMethod(shift=shift); cutoff=cutoff)
ω_krylov, _ = solve(solver_krylov, k_test; bands=1:4)
println("  KrylovKit: ω = ", round.(ω_krylov, digits=4))

# LOBPCG (similar to KrylovKit)
solver_lobpcg = Solver(FullVectorEM(), geo, resolution, LOBPCGMethod(shift=shift); cutoff=cutoff)
ω_lobpcg, _ = solve(solver_lobpcg, k_test; bands=1:4)
println("  LOBPCG:    ω = ", round.(ω_lobpcg, digits=4))

# Dense (for comparison, if system is small enough)
# Note: Dense computes ALL eigenvalues, filters ω² < shift, then sorts ascending.
#       This may give different bands than iterative methods which use shift-and-invert
#       to find eigenvalues closest to σ. Both results are valid eigenvalues.
if solver.basis.num_pw * 3 < 2000
    solver_dense = Solver(FullVectorEM(), geo, resolution, DenseMethod(shift=shift); cutoff=cutoff)
    ω_dense, _ = solve(solver_dense, k_test; bands=1:8)  # More bands to show all eigenvalues
    println("  Dense:     ω = ", round.(ω_dense, digits=4))
    println("  (Dense finds all eigenvalues in ascending order)")
end

# Plot band structure (scatter recommended for 3D due to Γ point anomalies)
println("\n=== Generating Plot ===")
p = plot_bands(bands_result;
    title="3D FCC Photonic Crystal (ε=12 spheres, r=0.25a)",
    ylabel="Frequency ω (2πc/a)",
    scatter=true,
    markersize=2,
    color=:blue,
    size=(700, 500)
)

# Save figure
output_file = joinpath(@__DIR__, "401_fcc_bands.png")
savefig(p, output_file)
println("Saved: $output_file")

println("\n✅ 3D FCC photonic crystal example complete!")
