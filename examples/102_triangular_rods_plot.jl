# Last-Modified: 2025-12-15T19:00:00+09:00
# Plot Triangular Lattice Photonic Crystal Band Structure
# Requires: Plots.jl

using PhoXonic
using Plots

println("Computing band structure...")

# Triangular (hexagonal) lattice
a = 1.0
lat = hexagonal_lattice(a)

# Geometry: air background with dielectric rod at center
air = Dielectric(1.0)
rod = Dielectric(12.0)
geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])

# Create solver for TM polarization
solver = Solver(TMWave(), geo, (64, 64); cutoff=7)

# Create k-path: Γ → M → K → Γ
kpath = simple_kpath_hexagonal(a=a, npoints=50)

# Compute bands
bands_result = compute_bands(solver, kpath; bands=1:8)

# Extract data
dists = bands_result.distances
freqs = bands_result.frequencies

# Create plot
p = plot(
    xlabel="Wave vector",
    ylabel="Frequency (ωa/2πc)",
    title="TM Bands: Triangular Lattice (ε=12 rods in air, r=0.2a)",
    legend=false,
    grid=true,
    size=(800, 500)
)

# Plot each band
for b in 1:size(freqs, 2)
    plot!(p, dists, freqs[:, b], linewidth=2, color=:blue)
end

# Add high-symmetry point labels
label_positions = [d for (i, _) in bands_result.labels for d in [dists[i]]]
label_names = [l for (_, l) in bands_result.labels]

vline!(p, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p, label_positions, label_names)

# Highlight band gap region (between band 1 and 2)
gap_info = find_bandgap(bands_result, 1, 2)
if gap_info.gap > 0
    hspan!(p, [gap_info.max_lower, gap_info.min_upper],
           alpha=0.2, color=:yellow, label="")
end

# Save plot
savefig(p, joinpath(@__DIR__, "102_triangular_tm_bands.png"))
println("Saved: 102_triangular_tm_bands.png")

# Show gap info
println("\nBand gap (bands 1-2): $(round(gap_info.gap_ratio*100, digits=1))% gap-to-midgap ratio")

display(p)
