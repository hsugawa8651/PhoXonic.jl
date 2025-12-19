# Last-Modified: 2025-12-18T20:00:00+09:00
# Simple example: Dense vs LOBPCG solver comparison
# Shows basic usage of different solvers

using PhoXonic
using Plots

println("=== Simple Solver Comparison ===")
println()

# ============================================================================
# Setup: Vasseur2001 Steel/Epoxy triangular lattice
# ============================================================================
ρ_steel = 7780.0
μ_steel = 8.1e10
λ_steel = 26.4e10 - 2*μ_steel

ρ_epoxy = 1142.0
μ_epoxy = 0.148e10
λ_epoxy = 0.754e10 - 2*μ_epoxy

steel = IsotropicElastic(ρ=ρ_steel, λ=λ_steel, μ=μ_steel)
epoxy = IsotropicElastic(ρ=ρ_epoxy, λ=λ_epoxy, μ=μ_epoxy)

a = 1.0
lat = hexagonal_lattice(a)
f_target = 0.4
r = sqrt(f_target * a^2 * sqrt(3) / (2π))

geo = Geometry(lat, epoxy, [(Circle([0.0, 0.0], r), steel)])

# ============================================================================
# Dense solver (default)
# ============================================================================
println("--- Dense Solver ---")
solver_dense = Solver(PSVWave(), geo, (64, 64); cutoff=10)
println("Plane waves: ", solver_dense.basis.num_pw)

kpath = simple_kpath_hexagonal(a=a, npoints=20)

t_dense = @elapsed begin
    bands_dense = compute_bands(solver_dense, kpath; bands=1:10, verbose=false)
end
println("Time: $(round(t_dense, digits=2)) s")

# ============================================================================
# LOBPCG solver (faster for large problems)
# ============================================================================
println("\n--- LOBPCG Solver ---")
solver_lobpcg = Solver(PSVWave(), geo, (64, 64), LOBPCGMethod(); cutoff=10)
println("Plane waves: ", solver_lobpcg.basis.num_pw)

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
println("Speedup: $(round(t_dense/t_lobpcg, digits=2))x")

# Check accuracy
max_diff = maximum(abs.(bands_dense.frequencies - bands_lobpcg.frequencies))
println("Max frequency difference: $(round(max_diff, digits=1)) rad/s")

# ============================================================================
# Plot comparison
# ============================================================================
p = plot(layout=(1,2), size=(900, 400), title=["Dense" "LOBPCG"])

for b in 1:size(bands_dense.frequencies, 2)
    plot!(p[1], bands_dense.distances, bands_dense.frequencies[:, b],
          color=:blue, legend=false)
    plot!(p[2], bands_lobpcg.distances, bands_lobpcg.frequencies[:, b],
          color=:red, legend=false)
end

xlabel!(p[1], "Wave Vector")
xlabel!(p[2], "Wave Vector")
ylabel!(p[1], "ω (rad/s)")

savefig(p, joinpath(@__DIR__, "207_solver_simple.png"))
println("\nSaved: 207_solver_simple.png")
