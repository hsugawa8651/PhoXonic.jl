# 210_vacuum_phononic.jl
# Phononic crystal with vacuum holes using ElasticVoid (Tanaka limit)
#
# Reference: Maldovan & Thomas, Appl. Phys. B 83, 595 (2006)
# DOI: https://doi.org/10.1007/s00340-006-2241-y
# Structure: Air holes in Si, triangular lattice, r/a = 0.46

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhoXonic
using Plots

println("=" ^ 60)
println("Phononic Crystal with Vacuum Holes (ElasticVoid)")
println("=" ^ 60)

# =============================================================================
# Material Parameters: Silicon with vacuum holes
# =============================================================================

# Silicon (from Maldovan 2006)
ρ_si = 2330.0   # kg/m³
c_T = 5360.0    # m/s (transverse velocity)
c_L = 8950.0    # m/s (longitudinal velocity)

# Convert to Lamé parameters
μ_si = ρ_si * c_T^2
λ_si = ρ_si * (c_L^2 - 2 * c_T^2)

si = IsotropicElastic(; ρ=ρ_si, λ=λ_si, μ=μ_si)

# Vacuum: ElasticVoid uses Tanaka limit (ρ/C → 0)
# Default ρ_ratio=1e-7 pushes spurious flat bands to high frequency
void = ElasticVoid()

println("\nMaterials:")
println("  Si: ρ=$(ρ_si) kg/m³, c_T=$(c_T) m/s, c_L=$(c_L) m/s")
println(
    "  Void: ρ=$(PhoXonic.density(void)), v_T=$(round(transverse_velocity(void), digits=1))"
)

# =============================================================================
# Geometry: Triangular lattice, r/a = 0.46
# =============================================================================

a = 1.0
r_over_a = 0.46

lat = hexagonal_lattice(a)
geo = Geometry(lat, si, [(Circle([0.0, 0.0], r_over_a * a), void)])

println("\nGeometry:")
println("  Lattice: Triangular (hexagonal)")
println("  r/a = $r_over_a")

# =============================================================================
# Compute band structures
# =============================================================================

resolution = (64, 64)
cutoff = 7
n_bands = 10

println("\nSolver: resolution=$resolution, cutoff=$cutoff")

kpath = simple_kpath_hexagonal(; npoints=50)

println("\nComputing bands...")
println("  SH...")
solver_sh = Solver(SHWave(), geo, resolution; cutoff=cutoff)
@time bands_sh = compute_bands(solver_sh, kpath; bands=1:n_bands)

println("  PSV...")
solver_psv = Solver(PSVWave(), geo, resolution; cutoff=cutoff)
@time bands_psv = compute_bands(solver_psv, kpath; bands=1:n_bands)

# =============================================================================
# Normalize and analyze
# =============================================================================

norm = 1.0 / (2π * c_T)
freqs_sh = bands_sh.frequencies * norm
freqs_psv = bands_psv.frequencies * norm

println("\n" * "=" ^ 60)
println("Band Gap Analysis (ωa/2πcT)")
println("=" ^ 60)

println("\nSH gaps:")
for gap in find_all_gaps(bands_sh)
    lower = gap.max_lower * norm
    upper = gap.min_upper * norm
    if lower < 2.0
        println(
            "  Band $(gap.bands): $(round(lower, digits=3)) - $(round(upper, digits=3))"
        )
    end
end

println("\nPSV gaps:")
for gap in find_all_gaps(bands_psv)
    lower = gap.max_lower * norm
    upper = gap.min_upper * norm
    if lower < 2.0
        println(
            "  Band $(gap.bands): $(round(lower, digits=3)) - $(round(upper, digits=3))"
        )
    end
end

println("\nExpected from Maldovan 2006: 0.830 - 0.963")

# =============================================================================
# Plot (SH and PSV as subplots)
# =============================================================================

dists = bands_sh.distances
labels = bands_sh.labels
label_pos = [dists[i] for (i, _) in labels]
label_names = [l for (_, l) in labels]

# SH subplot
p1 = plot(; ylabel="ωa/2πcT", title="SH Bands", legend=:topright, grid=true, ylims=(0, 1.5))
for b in 1:size(freqs_sh, 2)
    plot!(p1, dists, freqs_sh[:, b]; lw=2, color=:green, label="")
end
hspan!(p1, [0.830, 0.963]; alpha=0.2, color=:blue, label="Paper (complete)")
vline!(p1, label_pos; color=:gray, ls=:dash, alpha=0.5, label="")
xticks!(p1, label_pos, label_names)

# PSV subplot
p2 = plot(; title="PSV Bands", legend=:topright, grid=true, ylims=(0, 1.5))
for b in 1:size(freqs_psv, 2)
    plot!(p2, dists, freqs_psv[:, b]; lw=2, color=:orange, label="")
end
hspan!(p2, [0.830, 0.963]; alpha=0.2, color=:blue, label="Paper (complete)")
vline!(p2, label_pos; color=:gray, ls=:dash, alpha=0.5, label="")
xticks!(p2, label_pos, label_names)

# Combined
p = plot(
    p1,
    p2;
    layout=(1, 2),
    size=(1100, 450),
    plot_title="Maldovan 2006: Vacuum/Si (r/a=0.46)",
    left_margin=10Plots.mm,
)

savefig(p, joinpath(@__DIR__, "210_vacuum_phononic.png"))
println("\nSaved: 210_vacuum_phononic.png")
