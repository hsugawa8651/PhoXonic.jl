# Last-Modified: 2025-12-25T15:00:00+09:00
# Simple example: Dense vs LOBPCG solver comparison
#
# LOBPCG becomes faster than Dense for large problems (cutoff ≥ 12).
# With warm_start=true (default), LOBPCG reuses eigenvectors from
# previous k-points for faster convergence.

using PhoXonic
using Plots
default(
    guidefontsize=14,
    tickfontsize=12,
    titlefontsize=14,
    left_margin=10Plots.mm,
    right_margin=10Plots.mm,
    top_margin=5Plots.mm,
    bottom_margin=10Plots.mm,
)
using Printf

println("=== Simple Solver Comparison ===")
println()

# ============================================================================
# Setup: Vasseur2001 Steel/Epoxy triangular lattice
# ============================================================================
ρ_steel = 7780.0
μ_steel = 8.1e10
λ_steel = 26.4e10 - 2 * μ_steel

ρ_epoxy = 1142.0
μ_epoxy = 0.148e10
λ_epoxy = 0.754e10 - 2 * μ_epoxy

steel = IsotropicElastic(; ρ=ρ_steel, λ=λ_steel, μ=μ_steel)
epoxy = IsotropicElastic(; ρ=ρ_epoxy, λ=λ_epoxy, μ=μ_epoxy)

a = 1.0
lat = hexagonal_lattice(a)
f_target = 0.4
r = sqrt(f_target * a^2 * sqrt(3) / (2π))

geo = Geometry(lat, epoxy, [(Circle([0.0, 0.0], r), steel)])

# Use cutoff=15 where LOBPCG shows clear advantage
cutoff = 15
kpath = simple_kpath_hexagonal(; a=a, npoints=20)

# ============================================================================
# Dense solver (default, exact)
# ============================================================================
println("--- Dense Solver (reference) ---")
solver_dense = Solver(PSVWave(), geo, (64, 64); cutoff=cutoff)
n_pw = solver_dense.basis.num_pw
println("Plane waves: $n_pw (matrix size: $(2*n_pw)×$(2*n_pw))")

t_dense = @elapsed begin
    bands_dense = compute_bands(solver_dense, kpath; bands=1:10, verbose=false)
end
println("Time: $(round(t_dense, digits=2)) s")

# ============================================================================
# LOBPCG solver with warm start (default)
# ============================================================================
println("\n--- LOBPCG Solver (warm_start=true) ---")
solver_lobpcg = Solver(PSVWave(), geo, (64, 64), LOBPCGMethod(); cutoff=cutoff)
println("Plane waves: $n_pw (matrix size: $(2*n_pw)×$(2*n_pw))")

t_lobpcg = @elapsed begin
    bands_lobpcg = compute_bands(solver_lobpcg, kpath; bands=1:10, verbose=false)
end
println("Time: $(round(t_lobpcg, digits=2)) s")

# ============================================================================
# Results
# ============================================================================
println("\n=== Results ===")
println("Dense:  $(round(t_dense, digits=2)) s")
println("LOBPCG: $(round(t_lobpcg, digits=2)) s")
speedup = t_dense / t_lobpcg
println("Speedup: $(round(speedup, digits=2))x")

# Check accuracy
max_diff = maximum(abs.(bands_dense.frequencies - bands_lobpcg.frequencies))
freq_ref = max.(bands_dense.frequencies, 100.0)
max_rel_err = maximum(abs.(bands_dense.frequencies - bands_lobpcg.frequencies) ./ freq_ref)
println("Max frequency difference: $(round(max_diff, digits=1)) rad/s")
@printf("Max relative error: %.2f%%\n", max_rel_err * 100)

# ============================================================================
# Plot comparison
# ============================================================================
# Convert to krad/s for cleaner axis labels
freq_dense_k = bands_dense.frequencies / 1000
freq_lobpcg_k = bands_lobpcg.frequencies / 1000

ymax = max(maximum(freq_dense_k), maximum(freq_lobpcg_k)) * 1.05
ylims_common = (0, ymax)

# Get k-path labels
dists = bands_dense.distances
labels_data = bands_dense.labels
label_positions = [dists[i] for (i, _) in labels_data]
label_names = [l for (_, l) in labels_data]

title_dense = "Dense (cutoff=$cutoff, N=$(2*n_pw))"
title_lobpcg = @sprintf(
    "LOBPCG (%.1fx speedup)\nmax error: %.0f rad/s (%.2f%%)",
    speedup,
    max_diff,
    max_rel_err * 100
)

p = plot(; layout=(1, 2), size=(1200, 500), left_margin=5Plots.mm)

for b in 1:size(freq_dense_k, 2)
    plot!(
        p[1],
        dists,
        freq_dense_k[:, b];
        color=:blue,
        legend=false,
        ylims=ylims_common,
    )
    plot!(
        p[2],
        dists,
        freq_lobpcg_k[:, b];
        color=:red,
        legend=false,
        ylims=ylims_common,
    )
end

# Add k-path labels and vertical lines
for sp in [p[1], p[2]]
    vline!(sp, label_positions; color=:gray, linestyle=:dash, alpha=0.5)
    xticks!(sp, label_positions, label_names)
end

title!(p[1], title_dense)
title!(p[2], title_lobpcg)
xlabel!(p[1], "Wave Vector")
xlabel!(p[2], "Wave Vector")
ylabel!(p[1], "ω (10³ rad/s)")

savefig(p, joinpath(@__DIR__, "207_solver_simple.png"))
println("\nSaved: 207_solver_simple.png")
