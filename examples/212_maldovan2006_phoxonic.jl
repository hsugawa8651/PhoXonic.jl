# Last-Modified: 2025-12-21
# Maldovan & Thomas 2006: Si with Air Holes - Phoxonic Crystal
# Reference: M. Maldovan & E.L. Thomas, Appl. Phys. Lett. 88, 251907 (2006)
#
# Demonstrates:
# - Simultaneous photonic and phononic band structure calculation
# - ElasticVoid for vacuum regions in phononic crystals
# - Combined structure + dispersion figure for JOSS paper
#
# Structure: Triangular lattice, air holes in silicon, r/a = 0.46

using PhoXonic
using Printf
using Plots

println("=" ^ 60)
println("Maldovan & Thomas 2006: Phoxonic Crystal")
println("Si with Air Holes (r/a = 0.46)")
println("=" ^ 60)
println()

# ============================================================================
# Material properties
# ============================================================================

# Silicon
ε_si = 13.0              # permittivity
ρ_si = 2330.0            # kg/m³
c_T = 5360.0             # m/s (transverse velocity)
c_L = 8950.0             # m/s (longitudinal velocity)

# Lame parameters from velocities
μ_si = ρ_si * c_T^2
λ_si = ρ_si * (c_L^2 - 2 * c_T^2)

println("Silicon properties:")
println("  ε = $ε_si")
println("  ρ = $ρ_si kg/m³")
println("  c_T = $c_T m/s")
println("  c_L = $c_L m/s")

# ============================================================================
# Geometry: Triangular lattice with air holes in silicon
# ============================================================================

a = 1.0           # lattice constant (normalized)
r_over_a = 0.46   # radius ratio from paper

lat = hexagonal_lattice(a)

# Photonic: Si background with air holes
si_dielectric = Dielectric(ε_si)
air_dielectric = Dielectric(1.0)
geo_photonic = Geometry(lat, si_dielectric,
    [(Circle([0.0, 0.0], r_over_a * a), air_dielectric)])

# Phononic: Si background with void holes (ElasticVoid)
si_elastic = IsotropicElastic(ρ=ρ_si, λ=λ_si, μ=μ_si)
geo_phononic = Geometry(lat, si_elastic,
    [(Circle([0.0, 0.0], r_over_a * a), ElasticVoid())])

println("\nGeometry:")
println("  Lattice: Triangular (hexagonal)")
println("  Structure: Air holes in Silicon")
println("  r/a = $r_over_a")

# ============================================================================
# K-path: Γ → M → K → Γ
# ============================================================================

npoints = 40
kpath = simple_kpath_hexagonal(a=a, npoints=npoints)
println("\nK-path: Γ → M → K → Γ ($(length(kpath.points)) points)")

# ============================================================================
# Solver parameters
# ============================================================================

resolution = (64, 64)
cutoff = 7
n_bands_photonic = 6
n_bands_phononic = 8

println("\nSolver parameters:")
println("  Resolution: $resolution")
println("  Cutoff: $cutoff")

# ============================================================================
# Compute band structures
# ============================================================================

println("\n" * "=" ^ 60)
println("Computing band structures...")
println("=" ^ 60)

# Photonic
println("\nPhotonic:")
solver_te = Solver(TEWave(), geo_photonic, resolution; cutoff=cutoff)
solver_tm = Solver(TMWave(), geo_photonic, resolution; cutoff=cutoff)

print("  TE...")
@time bands_te = compute_bands(solver_te, kpath; bands=1:n_bands_photonic)
print("  TM...")
@time bands_tm = compute_bands(solver_tm, kpath; bands=1:n_bands_photonic)

# Phononic
println("\nPhononic:")
solver_sh = Solver(SHWave(), geo_phononic, resolution; cutoff=cutoff)
solver_psv = Solver(PSVWave(), geo_phononic, resolution; cutoff=cutoff)

print("  SH...")
@time bands_sh = compute_bands(solver_sh, kpath; bands=1:n_bands_phononic)
print("  PSV...")
@time bands_psv = compute_bands(solver_psv, kpath; bands=1:n_bands_phononic)

# ============================================================================
# Normalize frequencies
# ============================================================================

# Photonic: ωa/(2πc) - divide by √ε_background
photonic_norm = 1.0 / sqrt(ε_si)
freqs_te = bands_te.frequencies * photonic_norm
freqs_tm = bands_tm.frequencies * photonic_norm

# Phononic: ωa/(2πc_T)
phononic_norm = 1.0 / (2π * c_T)
freqs_sh = bands_sh.frequencies * phononic_norm
freqs_psv = bands_psv.frequencies * phononic_norm

# ============================================================================
# Find band gaps
# ============================================================================

println("\n" * "=" ^ 60)
println("Band Gap Analysis")
println("=" ^ 60)

gaps_te = find_all_gaps(bands_te; threshold=0.01)
gaps_tm = find_all_gaps(bands_tm; threshold=0.01)
gaps_sh = find_all_gaps(bands_sh; threshold=0.01)
gaps_psv = find_all_gaps(bands_psv; threshold=0.01)

println("\nTE gaps:")
for g in gaps_te
    println("  Bands $(g.bands): $(round(g.max_lower * photonic_norm, digits=3)) - $(round(g.min_upper * photonic_norm, digits=3))")
end

println("\nTM gaps:")
for g in gaps_tm
    println("  Bands $(g.bands): $(round(g.max_lower * photonic_norm, digits=3)) - $(round(g.min_upper * photonic_norm, digits=3))")
end

println("\nSH gaps:")
for g in gaps_sh
    println("  Bands $(g.bands): $(round(g.max_lower * phononic_norm, digits=3)) - $(round(g.min_upper * phononic_norm, digits=3))")
end

println("\nPSV gaps:")
for g in gaps_psv
    println("  Bands $(g.bands): $(round(g.max_lower * phononic_norm, digits=3)) - $(round(g.min_upper * phononic_norm, digits=3))")
end

# ============================================================================
# Create combined figure for JOSS paper
# ============================================================================

println("\n" * "=" ^ 60)
println("Creating JOSS paper figure...")
println("=" ^ 60)

dists = bands_te.distances
label_positions = [dists[i] for (i, _) in bands_te.labels]
label_names = [l for (_, l) in bands_te.labels]

# --- Structure plot ---
p_struct = plot(
    xlabel = "x/a",
    ylabel = "y/a",
    title = "Structure",
    aspect_ratio = :equal,
    xlims = (-0.6, 1.6),
    ylims = (-0.5, 1.1),
    legend = false,
    grid = false
)

# Draw Si background (extended unit cell region)
# Hexagonal unit cell vertices
uc_x = [0, 1, 1.5, 1, 0, -0.5, 0]
uc_y = [0, 0, sqrt(3)/2, sqrt(3), sqrt(3), sqrt(3)/2, 0]
plot!(p_struct, uc_x, uc_y,
    seriestype=:shape, fillcolor=:steelblue, linecolor=:black, linewidth=2,
    fillalpha=0.8, label="")

# Draw air holes (white circles)
θ = range(0, 2π, length=100)
r = r_over_a * a

# Hole positions for visualization
hole_positions = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.5, sqrt(3)/2],
    [-0.5, sqrt(3)/2],
    [1.5, sqrt(3)/2],
    [0.0, sqrt(3)],
    [1.0, sqrt(3)]
]

for pos in hole_positions
    cx, cy = pos
    hx = cx .+ r * cos.(θ)
    hy = cy .+ r * sin.(θ)
    plot!(p_struct, hx, hy, seriestype=:shape,
        fillcolor=:white, linecolor=:white, linewidth=0, label="")
end

# Unit cell outline (dashed, orange for visibility)
plot!(p_struct, uc_x, uc_y,
    linecolor=:darkorange, linewidth=2, linestyle=:dash, label="")

# Add "Air" label
annotate!(p_struct, 0.0, 0.0, text("Air", :black, 10))

# --- TE plot ---
p_te = plot(
    xlabel = "",
    ylabel = "ωa/2πc",
    title = "TE (photonic)",
    legend = false,
    grid = true,
    ylims = (0, 0.8),
    xticks = (label_positions, label_names)
)

for b in 1:size(freqs_te, 2)
    plot!(p_te, dists, freqs_te[:, b], linewidth=2, color=:blue)
end

# TE gap shading
for g in gaps_te
    gmin = g.max_lower * photonic_norm
    gmax = g.min_upper * photonic_norm
    hspan!(p_te, [gmin, gmax], alpha=0.2, color=:green, label="")
end

vline!(p_te, label_positions, color=:gray, linestyle=:dash, alpha=0.5)

# --- TM plot ---
p_tm = plot(
    xlabel = "",
    ylabel = "ωa/2πc",
    title = "TM (photonic)",
    legend = false,
    grid = true,
    ylims = (0, 0.8),
    xticks = (label_positions, label_names)
)

for b in 1:size(freqs_tm, 2)
    plot!(p_tm, dists, freqs_tm[:, b], linewidth=2, color=:red)
end

# TM gap shading
for g in gaps_tm
    gmin = g.max_lower * photonic_norm
    gmax = g.min_upper * photonic_norm
    hspan!(p_tm, [gmin, gmax], alpha=0.2, color=:green, label="")
end

vline!(p_tm, label_positions, color=:gray, linestyle=:dash, alpha=0.5)

# --- SH plot ---
p_sh = plot(
    xlabel = "Wave vector",
    ylabel = "ωa/2πcT",
    title = "SH (phononic)",
    legend = false,
    grid = true,
    ylims = (0, 1.5),
    xticks = (label_positions, label_names)
)

for b in 1:size(freqs_sh, 2)
    plot!(p_sh, dists, freqs_sh[:, b], linewidth=2, color=:green)
end

# SH gap shading
for g in gaps_sh
    gmin = g.max_lower * phononic_norm
    gmax = g.min_upper * phononic_norm
    hspan!(p_sh, [gmin, gmax], alpha=0.2, color=:orange, label="")
end

vline!(p_sh, label_positions, color=:gray, linestyle=:dash, alpha=0.5)

# --- PSV plot ---
p_psv = plot(
    xlabel = "Wave vector",
    ylabel = "ωa/2πcT",
    title = "PSV (phononic)",
    legend = false,
    grid = true,
    ylims = (0, 1.5),
    xticks = (label_positions, label_names)
)

for b in 1:size(freqs_psv, 2)
    plot!(p_psv, dists, freqs_psv[:, b], linewidth=2, color=:orange)
end

# PSV gap shading
for g in gaps_psv
    gmin = g.max_lower * phononic_norm
    gmax = g.min_upper * phononic_norm
    hspan!(p_psv, [gmin, gmax], alpha=0.2, color=:orange, label="")
end

vline!(p_psv, label_positions, color=:gray, linestyle=:dash, alpha=0.5)

# --- Combined layout ---
# Layout: [structure] [TE] [TM]
#                     [SH] [PSV]

l = @layout [
    a{0.35w} [b c
              d e]
]

p_combined = plot(p_struct, p_te, p_tm, p_sh, p_psv,
    layout = l,
    size = (1400, 700),
    plot_title = "Maldovan 2006: PhoXonic Crystal - Si with Air Holes (r/a=0.46)",
    plot_titlefontsize = 14
)

# Save
savefig(p_combined, joinpath(@__DIR__, "212_maldovan2006_phoxonic.png"))
println("\nSaved: 212_maldovan2006_phoxonic.png")

display(p_combined)

println("\n" * "=" ^ 60)
println("Expected results from paper (r/a = 0.46):")
println("  Photonic gap: ωa/2πc ≈ 0.408 - 0.464")
println("  Phononic gap: ωa/2πcT ≈ 0.830 - 0.963")
println("=" ^ 60)
println("\nDone!")
