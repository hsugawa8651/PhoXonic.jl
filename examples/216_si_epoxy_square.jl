# Last-Modified: 2026-02-24
# Si/Epoxy phononic crystal band structure (square lattice)
# Reference: Dobrzynski et al., "Phononics", Elsevier (2017), Ch.5, Fig. 20A
# ISBN: 978-0-12-809931-5
# Material parameters from Table 2, p.300
#
# Structure: Circular Si cylinders in Epoxy matrix, square lattice
# Filling fraction: f = 0.68

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
println("Dobrzynski 2017, Ch.5 Fig. 20A: Si/Epoxy Phononic Crystal")
println("Square lattice, Si cylinders in Epoxy matrix")
println("=" ^ 60)
println()

# ============================================================================
# Material properties (Table 2, p.300)
# ============================================================================
# Silicon (inclusion) - hard material
# C11 = 165.7 GPa, C44 = 79.62 GPa, C12 = 63.9 GPa
ρ_si = 2331.0         # kg/m³
λ_si = 6.46e9         # Pa (C12 = C11 - 2*C44 → λ = C12 = 63.9? No: λ = C11-2μ = 165.7-2*79.62 = 6.46 GPa)
μ_si = 79.62e9        # Pa (C44)

# Epoxy (matrix) - soft material
# C11 = 7.61 GPa, C44 = 1.59 GPa, C12 = 4.43 GPa
ρ_epoxy = 1180.0      # kg/m³
λ_epoxy = 4.43e9      # Pa (C11 - 2*C44 = 7.61 - 2*1.59 = 4.43 GPa)
μ_epoxy = 1.59e9      # Pa (C44)

silicon = IsotropicElastic(; ρ=ρ_si, λ=λ_si, μ=μ_si)
epoxy = IsotropicElastic(; ρ=ρ_epoxy, λ=λ_epoxy, μ=μ_epoxy)

println("Materials:")
println(
    "  Si:    ρ = $ρ_si kg/m³, C44 = $(silicon.C44/1e9) GPa, cₜ = $(round(sqrt(silicon.C44/silicon.ρ), digits=0)) m/s",
)
println(
    "  Epoxy: ρ = $ρ_epoxy kg/m³, C44 = $(epoxy.C44/1e9) GPa, cₜ = $(round(sqrt(epoxy.C44/epoxy.ρ), digits=0)) m/s",
)

# ============================================================================
# Geometry: square lattice with Si cylinders
# ============================================================================
a = 1.0  # lattice constant (normalized)
lat = square_lattice(a)

f = 0.68  # filling fraction
r = a * sqrt(f / π)  # r ≈ 0.465a
geo = Geometry(lat, epoxy, [(Circle([0.0, 0.0], r), silicon)])

println("\nGeometry:")
println("  Square lattice, a = $a")
println("  Cylinder radius r = $(round(r, digits=4))")
println("  Filling fraction f = $(round(π*r^2/a^2, digits=2))")

# ============================================================================
# K-path: M → Γ → X → M (matching book's Fig. 20A)
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
nbands_sh = 1:8
nbands_psv = 1:12

# ============================================================================
# SH Wave
# ============================================================================
println("\nComputing SH bands...")
solver_sh = Solver(SHWave(), geo, (64, 64); cutoff=cutoff)
t_sh = @elapsed bands_sh = compute_bands(solver_sh, kpath; bands=nbands_sh)
println("  Done in $(round(t_sh, digits=2)) s")

# ============================================================================
# P-SV Wave
# ============================================================================
println("Computing P-SV bands...")
solver_psv = Solver(PSVWave(), geo, (64, 64); cutoff=cutoff)
t_psv = @elapsed bands_psv = compute_bands(solver_psv, kpath; bands=nbands_psv)
println("  Done in $(round(t_psv, digits=2)) s")

# ============================================================================
# Plot band structure
# ============================================================================
println("\nCreating plots...")

# Normalize to fa (m/s): fa = ω/(2π) × a, matching book's y-axis
fa_factor = a / (2π)

freqs_sh_fa = bands_sh.frequencies * fa_factor
freqs_psv_fa = bands_psv.frequencies * fa_factor

dists = bands_sh.distances
label_positions = [dists[i] for (i, _) in bands_sh.labels]
label_names = [l for (_, l) in bands_sh.labels]

# Combined: SH (dashed) + PSV (solid) on same plot, like Fig. 20A
p = plot(;
    xlabel="Wave vector",
    ylabel="fa (m/s)",
    title="Dobrzynski 2017 Ch.5 Fig.20A: Si/Epoxy Square (f=0.68)",
    legend=false,
    grid=true,
    ylim=(0, 4000),
    size=(700, 500),
)
for b in 1:size(freqs_sh_fa, 2)
    plot!(p, dists, freqs_sh_fa[:, b]; linewidth=2, color=:blue, linestyle=:dash)
end
for b in 1:size(freqs_psv_fa, 2)
    plot!(p, dists, freqs_psv_fa[:, b]; linewidth=2, color=:red)
end
vline!(p, label_positions; color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p, label_positions, label_names)

savefig(p, joinpath(@__DIR__, "216_si_epoxy_square.png"))
println("Saved: 216_si_epoxy_square.png")

display(p)
println("\nDone!")
