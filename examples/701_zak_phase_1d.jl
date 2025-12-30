# Example 701: 1D Zak Phase Calculation
#
# This example demonstrates calculating the Zak phase for a 1D photonic crystal.
# The Zak phase is quantized to 0 or π for systems with inversion symmetry.

using PhoXonic

# Create a 1D photonic crystal with a dielectric layer
# Unit cell: air-dielectric-air
lat = lattice_1d(1.0)

# Geometry A: dielectric layer centered in unit cell
geo_A = Geometry(lat, Dielectric(1.0), [(Segment(0.35, 0.65), Dielectric(9.0))])

# Geometry B: dielectric layer at edges (same material, different origin)
geo_B = Geometry(lat, Dielectric(1.0), [
    (Segment(0.0, 0.15), Dielectric(9.0)),
    (Segment(0.85, 1.0), Dielectric(9.0)),
])

# Create solvers
solver_A = Solver(Photonic1D(), geo_A, 128; cutoff = 20)
solver_B = Solver(Photonic1D(), geo_B, 128; cutoff = 20)

# Compute Zak phases for bands 1-3
println("Computing Zak phases...")
println()

zak_A = compute_zak_phase(solver_A, 1:3; n_k = 100)
zak_B = compute_zak_phase(solver_B, 1:3; n_k = 100)

println("Geometry A (centered dielectric):")
for (i, phase) in enumerate(zak_A.phases)
    phase_norm = phase / π
    println("  Band $i: Zak phase = $(round(phase, digits=4)) ($(round(phase_norm, digits=2))π)")
end
println()

println("Geometry B (edge dielectric):")
for (i, phase) in enumerate(zak_B.phases)
    phase_norm = phase / π
    println("  Band $i: Zak phase = $(round(phase, digits=4)) ($(round(phase_norm, digits=2))π)")
end
println()

# The Zak phase should be quantized to 0 or π due to inversion symmetry
# Different unit cell choices can give different Zak phases
println("Note: For inversion-symmetric systems, Zak phase is quantized to 0 or π.")
println("      Different unit cell choices yield different Zak phases.")
