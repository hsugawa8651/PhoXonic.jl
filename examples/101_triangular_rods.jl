# Last-Modified: 2025-12-15T19:00:00+09:00
# Triangular Lattice Photonic Crystal Band Structure
# Dielectric rods in air - TM polarization

using PhoXonic
using LinearAlgebra

println("=== Triangular Lattice Photonic Crystal ===")
println("Dielectric rods (ε=12) in air, r=0.2a")
println()

# Triangular (hexagonal) lattice
a = 1.0
lat = hexagonal_lattice(a)
println("Lattice vectors:")
println("  a1 = ", lat.vectors[1])
println("  a2 = ", lat.vectors[2])

# Geometry: air background with dielectric rod at center
air = Dielectric(1.0)
rod = Dielectric(12.0)  # High dielectric like GaAs
geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])

# Create solver for TM polarization (this has a gap for rods in air)
println("\nCreating TM solver (cutoff=7, resolution=64x64)...")
solver = Solver(TMWave(), geo, (64, 64); cutoff=7)
println("Number of plane waves: ", solver.basis.num_pw)

# Create k-path: Γ → M → K → Γ
println("\nComputing band structure along Γ → M → K → Γ...")
kpath = simple_kpath_hexagonal(; a=a, npoints=30)
println("Number of k-points: ", length(kpath))

# Compute bands
bands_result = compute_bands(solver, kpath; bands=1:8, verbose=false)
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
            round(g.gap_ratio*100; digits=1),
            "% gap-to-midgap)",
        )
    end
end

# Output data for plotting
println("\n=== Band Structure Data ===")
println("# distance  band1  band2  band3  band4  band5  band6  band7  band8")
freqs = bands_result.frequencies
dists = bands_result.distances
for i in 1:length(dists)
    print(round(dists[i]; digits=4))
    for b in 1:8
        print("  ", round(freqs[i, b]; digits=4))
    end
    println()
end

# Print label positions
println("\n=== High Symmetry Points ===")
for (idx, label) in bands_result.labels
    println(label, " at distance ", round(dists[idx]; digits=4))
end
