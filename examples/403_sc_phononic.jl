# Last-Modified: 2025-12-15T19:00:00+09:00
# 3D Simple Cubic Phononic Crystal Band Structure
# Steel spheres in epoxy matrix - Full elastic waves

using PhoXonic
using LinearAlgebra
using Plots

println("=== 3D Simple Cubic Phononic Crystal ===")
println("Steel spheres in epoxy matrix")
println()

# Simple cubic lattice
a = 1.0  # lattice constant (normalized)
lat = cubic_lattice(a)
println("Simple cubic lattice (a = $a):")
println("  Lattice vectors:")
for (i, v) in enumerate(lat.vectors)
    println("    a$i = ", round.(v, digits=4))
end

# Materials (normalized units)
# Using Lamé parameters: λ (first Lamé), μ (shear modulus)
# Relations: C11 = λ + 2μ, C44 = μ

# Steel: high stiffness, high density
# Typical values (normalized): ρ=7.8, λ=120, μ=80 → C11=280
steel = IsotropicElastic(ρ=7.8, λ=120.0, μ=80.0)

# Epoxy: lower stiffness, lower density
# Typical values (normalized): ρ=1.2, λ=5, μ=1.5 → C11=8
epoxy = IsotropicElastic(ρ=1.2, λ=5.0, μ=1.5)

# Alternatively, create from Young's modulus and Poisson's ratio:
# steel = from_E_ν(7.8, 200.0, 0.3)
# epoxy = from_E_ν(1.2, 4.0, 0.35)

radius = 0.3  # sphere radius in units of a
geo = Geometry(lat, epoxy, [(Sphere([0.0, 0.0, 0.0], radius), steel)])

println("\nMaterials:")
println("  Steel sphere: ρ=$(steel.ρ), C11=$(steel.C11), C44=$(steel.C44)")
println("  Epoxy matrix: ρ=$(epoxy.ρ), C11=$(epoxy.C11), C44=$(epoxy.C44)")
println("  Sphere radius: r/a = $radius")

# Wave velocities for reference
v_L_steel = sqrt(steel.C11 / steel.ρ)
v_T_steel = sqrt(steel.C44 / steel.ρ)
v_L_epoxy = sqrt(epoxy.C11 / epoxy.ρ)
v_T_epoxy = sqrt(epoxy.C44 / epoxy.ρ)
println("\nWave velocities:")
println("  Steel: v_L = $(round(v_L_steel, digits=2)), v_T = $(round(v_T_steel, digits=2))")
println("  Epoxy: v_L = $(round(v_L_epoxy, digits=2)), v_T = $(round(v_T_epoxy, digits=2))")

# Create solver with shift-and-invert
# For 3D elastic, we need shift to target low-frequency acoustic modes
println("\nCreating FullElastic solver...")
resolution = (12, 12, 12)
cutoff = 7
shift = 0.1  # target low frequencies

solver = Solver(FullElastic(), geo, resolution, KrylovKitMethod(shift=shift); cutoff=cutoff)
println("  Resolution: $resolution")
println("  Cutoff: $cutoff")
println("  Number of plane waves: ", solver.basis.num_pw)
println("  Total DOF: ", solver.basis.num_pw * 3, " (3 displacement components)")
println("  Shift (σ): $shift")

# Create k-path: Γ → X → M → Γ → R → X (standard SC path)
println("\nComputing band structure...")
kpath = simple_kpath_cubic(a=a, npoints=15)
println("K-path: Γ → X → M → Γ → R → X")
println("Number of k-points: ", length(kpath))
println("Note: Γ point has ω=0 for acoustic modes (3 branches)")

# Compute bands
nbands = 8  # 3 acoustic + several optical branches
println("\nComputing $nbands bands...")
println("This may take a few minutes for 3D elastic calculations...")

bands_result = compute_bands(solver, kpath; bands=1:nbands, verbose=true)
println("Computation complete.")

# Find band gaps
println("\n=== Band Gaps ===")
gaps = find_all_gaps(bands_result; threshold=0.01)
if isempty(gaps)
    println("No significant band gaps found.")
else
    for g in gaps
        println("Gap between bands ", g.bands, ": Δω = ", round(g.gap, digits=4),
                " (", round(g.gap_ratio*100, digits=1), "% gap-to-midgap)")
    end
end

# Output key frequencies
println("\n=== Band Structure Summary ===")
freqs = bands_result.frequencies
dists = bands_result.distances

# Find X point (first high symmetry point after Γ)
x_idx = findfirst(idx -> idx > 1, [idx for (idx, _) in bands_result.labels])
if x_idx !== nothing
    x_point_idx = bands_result.labels[x_idx][1]
    println("At X point:")
    for b in 1:min(nbands, 6)
        println("  Band $b: ω = ", round(freqs[x_point_idx, b], digits=4))
    end
end

# Print label positions
println("\n=== High Symmetry Points ===")
for (idx, label) in bands_result.labels
    println(label, " at distance ", round(dists[idx], digits=4))
end

# Compare with different solver methods at a test point
println("\n=== Solver Method Comparison ===")
k_test = [0.5, 0.0, 0.0]  # X point
println("Testing at k = $k_test (X point)")

# KrylovKit
ω_krylov, _ = solve(solver, k_test; bands=1:4)
println("  KrylovKit: ω = ", round.(ω_krylov, digits=4))

# LOBPCG
solver_lobpcg = Solver(FullElastic(), geo, resolution, LOBPCGMethod(shift=shift); cutoff=cutoff)
ω_lobpcg, _ = solve(solver_lobpcg, k_test; bands=1:4)
println("  LOBPCG:    ω = ", round.(ω_lobpcg, digits=4))

# Plot band structure
println("\n=== Generating Plot ===")
p = plot_bands(bands_result;
    title="3D SC Phononic Crystal (Steel spheres in Epoxy, r=0.3a)",
    ylabel="Frequency ω (normalized)",
    scatter=true,
    markersize=2,
    color=:red,
    size=(700, 500)
)

# Save figure
output_file = joinpath(@__DIR__, "403_sc_phononic_bands.png")
savefig(p, output_file)
println("Saved: $output_file")

println("\n=== 3D phononic crystal example complete! ===")
