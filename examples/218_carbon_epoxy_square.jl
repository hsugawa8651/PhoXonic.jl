# Last-Modified: 2026-02-24
# Carbon/Epoxy phononic crystal band structure
# Reference: Dobrzynski et al., "Phononics", Elsevier (2017), Ch.5, Fig. 2A
# ISBN: 978-0-12-809931-5
# Material parameters from Table 1, p.276
# Normalization: Ω = ωa/(2πc₀), c₀ = √(C̄₄₄/ρ̄) (volume-averaged, p.277)
#
# Structure: Circular Carbon cylinders in Epoxy matrix, square lattice
# Filling fraction: f = 0.55

using PhoXonic
using Plots
default(;
    guidefontsize=14,
    tickfontsize=12,
    titlefontsize=14,
    left_margin=10Plots.mm,
    right_margin=10Plots.mm,
    top_margin=5Plots.mm,
    bottom_margin=10Plots.mm,
)

println("=" ^ 60)
println("Dobrzynski 2017, Ch.5 Fig. 2A: Carbon/Epoxy Phononic Crystal")
println("Square lattice, Carbon cylinders in Epoxy matrix")
println("=" ^ 60)
println()

# ============================================================================
# Material properties (Table 1, p.276)
# ============================================================================
# Carbon (inclusion)
ρ_carbon = 1750.0     # kg/m³
λ_carbon = 132.7e9    # Pa (C11 - 2*C44 = 309.6 - 2*88.46 = 132.7 GPa)
μ_carbon = 88.46e9    # Pa (C44)

# Epoxy (matrix)
ρ_epoxy = 1200.0      # kg/m³
λ_epoxy = 6.38e9      # Pa (C11 - 2*C44 = 9.64 - 2*1.63 = 6.38 GPa)
μ_epoxy = 1.63e9      # Pa (C44)

carbon = IsotropicElastic(; ρ=ρ_carbon, λ=λ_carbon, μ=μ_carbon)
epoxy = IsotropicElastic(; ρ=ρ_epoxy, λ=λ_epoxy, μ=μ_epoxy)

println("Materials:")
println("  Carbon: ρ = $ρ_carbon kg/m³, C44 = $(carbon.C44/1e9) GPa")
println("  Epoxy:  ρ = $ρ_epoxy kg/m³, C44 = $(epoxy.C44/1e9) GPa")

# ============================================================================
# Geometry: square lattice with carbon cylinders
# ============================================================================
a = 1.0  # lattice constant (normalized)
lat = square_lattice(a)

f = 0.55  # filling fraction
r = a * sqrt(f / π)  # r ≈ 0.418a
geo = Geometry(lat, epoxy, [(Circle([0.0, 0.0], r), carbon)])

println("\nGeometry:")
println("  Square lattice, a = $a")
println("  Cylinder radius r = $(round(r, digits=4))")
println("  Filling fraction f = $(round(π*r^2/a^2, digits=2))")

# ============================================================================
# K-path: M → Γ → X → M (matching book's Fig. 2A)
# ============================================================================
npoints = 40
scale = 2π / a
M = [0.5 * scale, 0.5 * scale]
Γ = [0.0, 0.0]
X = [0.5 * scale, 0.0]
kpath = SimpleKPath([M, Γ, X, M], ["M", "Γ", "X", "M"]; npoints=npoints)
println("\nK-path: $(length(kpath.points)) points")

# ============================================================================
# Solver parameters
# ============================================================================
cutoff = 10
nbands_sh = 1:15
nbands_psv = 1:15

println("\nSolver parameters:")
println("  Cutoff = $cutoff")
println("  SH bands = $nbands_sh")
println("  P-SV bands = $nbands_psv")

# ============================================================================
# SH Wave (scalar, out-of-plane = Z mode)
# ============================================================================
println("\n" * "=" ^ 60)
println("SH Wave (out-of-plane, Z mode)")
println("=" ^ 60)

solver_sh = Solver(SHWave(), geo, (64, 64); cutoff=cutoff)
println("Matrix size: $(matrix_dimension(solver_sh))")

t_sh = @elapsed bands_sh = compute_bands(solver_sh, kpath; bands=nbands_sh)
println("Dense: $(round(t_sh, digits=2)) s")

# ============================================================================
# P-SV Wave (vector, in-plane = XY mode)
# ============================================================================
println("\n" * "=" ^ 60)
println("P-SV Wave (in-plane, XY mode)")
println("=" ^ 60)

solver_psv = Solver(PSVWave(), geo, (64, 64); cutoff=cutoff)
println("Matrix size: $(matrix_dimension(solver_psv))")

t_psv = @elapsed bands_psv = compute_bands(solver_psv, kpath; bands=nbands_psv)
println("Dense: $(round(t_psv, digits=2)) s")

# ============================================================================
# Plot band structure
# ============================================================================
println("\nCreating plots...")

# Normalize: Ω = ωa/(2πc₀), c₀ = √(C̄₄₄/ρ̄) (volume-averaged, p.277)
ρ_bar = ρ_carbon * f + ρ_epoxy * (1 - f)
C44_bar = μ_carbon * f + μ_epoxy * (1 - f)
c0 = sqrt(C44_bar / ρ_bar)
norm_factor = a / (2π * c0)

println("  ρ̄ = $(round(ρ_bar, digits=1)) kg/m³")
println("  C̄₄₄ = $(round(C44_bar/1e9, digits=2)) GPa")
println("  c₀ = $(round(c0, digits=1)) m/s")

freqs_sh_norm = bands_sh.frequencies * norm_factor
freqs_psv_norm = bands_psv.frequencies * norm_factor

dists = bands_sh.distances
label_positions = [dists[i] for (i, _) in bands_sh.labels]
label_names = [l for (_, l) in bands_sh.labels]

# Combined plot: Z (dashed) + XY (solid) on same plot, like Fig. 2A
p = plot(;
    xlabel="Reduced wave vector",
    ylabel="Reduced frequency",
    title="Dobrzynski 2017 Ch.5 Fig.2A: C/Epoxy Square (f=0.55)",
    legend=false,
    grid=true,
    ylim=(0, 0.65),
    size=(700, 600),
)
# XY modes (solid) = P-SV
for b in 1:size(freqs_psv_norm, 2)
    plot!(p, dists, freqs_psv_norm[:, b]; linewidth=1.5, color=:black, linestyle=:solid)
end
# Z modes (dashed) = SH
for b in 1:size(freqs_sh_norm, 2)
    plot!(p, dists, freqs_sh_norm[:, b]; linewidth=1.5, color=:black, linestyle=:dash)
end
vline!(p, label_positions; color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p, label_positions, label_names)

savefig(p, joinpath(@__DIR__, "218_carbon_epoxy_square.png"))
println("\nSaved: 218_carbon_epoxy_square.png")

display(p)
println("\nDone!")
