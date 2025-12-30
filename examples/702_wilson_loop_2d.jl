# Example 702: 2D Wilson Loop Spectrum
#
# This example demonstrates calculating the Wilson loop spectrum for a 2D photonic crystal.
# The Wilson spectrum reveals topological properties through the winding of phases.

using PhoXonic

# Create a 2D photonic crystal: square lattice with circular rod
lat = square_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Circle([0.0, 0.0], 0.3), Dielectric(9.0))])

# Create solver for TM polarization
solver = Solver(TMWave(), geo, (32, 32); cutoff = 5)

println("Computing Wilson loop spectrum for 2D photonic crystal...")
println("Structure: square lattice, circular rod (r=0.3a), ε=9.0")
println()

# Compute Wilson spectrum for bands 1-2
# Loop direction :b2 means we scan along b1 and compute loops in b2 direction
result = compute_wilson_spectrum(solver, 1:2; n_k_path = 21, n_k_loop = 50, loop_direction = :b2)

println("Wilson spectrum calculated:")
println("  k-path points: $(length(result.k_values))")
println("  Number of bands: $(length(result.bands))")
println("  Loop direction: $(result.loop_direction)")
println()

# Display phases at a few k-points
println("Wilson phases at selected k-points:")
for i in [1, 6, 11, 16, 21]
    k = result.k_values[i]
    phases = result.phases[i, :]
    println("  k = $(round(k, digits=2)): phases = $(round.(phases, digits=3))")
end
println()

# Calculate winding numbers
println("Winding numbers:")
for band in 1:length(result.bands)
    w = winding_number(result, band)
    println("  Band $band: winding = $w")
end
println()

# For a simple photonic crystal, winding should be 0 (trivial topology)
println("Note: Non-zero winding indicates non-trivial topology.")
println("      This simple structure has trivial topology (winding = 0).")
