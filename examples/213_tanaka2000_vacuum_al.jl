# 213_tanaka2000_vacuum_al.jl
# Reproduce Tanaka et al. 2000, Phys. Rev. B 62, 7387
# DOI: https://doi.org/10.1103/PhysRevB.62.7387
# Figure 4: Vacuum cylinders in Al, square lattice, f=0.55
#
# This is the original paper that proposed the Tanaka limit (ρ/C → 0)
# for handling void regions in PWE calculations.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhoXonic
using Plots
default(guidefontsize=14, tickfontsize=12, titlefontsize=14, left_margin=10Plots.mm, right_margin=10Plots.mm, top_margin=5Plots.mm, bottom_margin=10Plots.mm)

println("=" ^ 60)
println("Tanaka et al. 2000 - Vacuum/Al Square Lattice")
println("=" ^ 60)

# =============================================================================
# Material Parameters from Tanaka 2000
# =============================================================================

# Aluminum
ρ_al = 2700.0    # kg/m³
c_T_al = 3110.0  # m/s (transverse)
c_L_al = 6420.0  # m/s (longitudinal)

# Convert to Lamé parameters
μ_al = ρ_al * c_T_al^2
λ_al = ρ_al * (c_L_al^2 - 2 * c_T_al^2)

al = IsotropicElastic(; ρ=ρ_al, λ=λ_al, μ=μ_al)
void = ElasticVoid()

println("\nMaterials:")
println("  Al: ρ=$(ρ_al) kg/m³, c_T=$(c_T_al) m/s, c_L=$(c_L_al) m/s")
println("  Void: ElasticVoid (Tanaka limit)")

# =============================================================================
# Geometry: Square lattice, filling fraction f=0.55
# =============================================================================

a = 1.0
f = 0.55  # filling fraction of vacuum
r = sqrt(f / π) * a  # radius from πr²/a² = f

lat = square_lattice(a)
geo = Geometry(lat, al, [(Circle([0.0, 0.0], r), void)])

println("\nGeometry:")
println("  Lattice: Square")
println("  Filling fraction f = $f")
println("  r/a = $(round(r/a, digits=3))")

# =============================================================================
# Compute band structures
# =============================================================================

resolution = (64, 64)
cutoff = 7
n_bands = 10

println("\nSolver: resolution=$resolution, cutoff=$cutoff")

kpath = simple_kpath_square(; npoints=50)

println("\nComputing bands...")
println("  SH...")
solver_sh = Solver(SHWave(), geo, resolution; cutoff=cutoff)
@time bands_sh = compute_bands(solver_sh, kpath; bands=1:n_bands)

println("  PSV...")
solver_psv = Solver(PSVWave(), geo, resolution; cutoff=cutoff)
@time bands_psv = compute_bands(solver_psv, kpath; bands=1:n_bands)

# =============================================================================
# Normalize: ωa/v_T (to match Tanaka 2000 Fig. 4)
# Note: Paper uses ωa/v_T, not ωa/2πc_T
# =============================================================================

norm = a / c_T_al  # No 2π to match paper
freqs_sh = bands_sh.frequencies * norm
freqs_psv = bands_psv.frequencies * norm

println("\n" * "=" ^ 60)
println("Band Gap Analysis (ωa/vT)")
println("=" ^ 60)

println("\nSH gaps:")
for gap in find_all_gaps(bands_sh)
    lower = gap.max_lower * norm
    upper = gap.min_upper * norm
    if lower < 6.0
        println(
            "  Band $(gap.bands): $(round(lower, digits=3)) - $(round(upper, digits=3)), ratio=$(round(gap.gap_ratio*100, digits=1))%",
        )
    end
end

println("\nPSV gaps:")
for gap in find_all_gaps(bands_psv)
    lower = gap.max_lower * norm
    upper = gap.min_upper * norm
    if lower < 6.0
        println(
            "  Band $(gap.bands): $(round(lower, digits=3)) - $(round(upper, digits=3)), ratio=$(round(gap.gap_ratio*100, digits=1))%",
        )
    end
end

# =============================================================================
# Plot (SH and PSV as subplots)
# =============================================================================

dists = bands_sh.distances
labels = bands_sh.labels
label_pos = [dists[i] for (i, _) in labels]
label_names = [l for (_, l) in labels]

# SH subplot
p1 = plot(; ylabel="ωa/vT", title="SH Bands", legend=false, grid=true, ylims=(0, 6))
for b in 1:size(freqs_sh, 2)
    plot!(p1, dists, freqs_sh[:, b]; lw=2, color=:green)
end
vline!(p1, label_pos; color=:gray, ls=:dash, alpha=0.5)
xticks!(p1, label_pos, label_names)

# PSV subplot
p2 = plot(; title="PSV Bands", legend=false, grid=true, ylims=(0, 6))
for b in 1:size(freqs_psv, 2)
    plot!(p2, dists, freqs_psv[:, b]; lw=2, color=:orange)
end
vline!(p2, label_pos; color=:gray, ls=:dash, alpha=0.5)
xticks!(p2, label_pos, label_names)

# Combined
p = plot(
    p1,
    p2;
    layout=(1, 2),
    size=(1200, 500),
    plot_title="Tanaka 2000 Fig.4: Vacuum/Al (f=0.55)",
    left_margin=10Plots.mm,
)

savefig(p, joinpath(@__DIR__, "213_tanaka2000_vacuum_al.png"))
println("\nSaved: 213_tanaka2000_vacuum_al.png")

println("\n" * "=" ^ 60)
println("Reference: Tanaka et al., Phys. Rev. B 62, 7387 (2000)")
println("=" ^ 60)
