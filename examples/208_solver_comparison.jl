# Last-Modified: 2025-12-18T20:00:00+09:00
# Detailed solver comparison: Dense vs LOBPCG
# Compares performance and accuracy across different problem sizes

using PhoXonic
using Plots
using Printf

println("=== Detailed Solver Comparison ===")
println("Dense vs LOBPCG performance analysis")
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

# k-path with moderate resolution
kpath = simple_kpath_hexagonal(a=a, npoints=30)

# ============================================================================
# Benchmark across different cutoff values
# ============================================================================
cutoffs = [5, 8, 10, 12, 15]
nbands = 15

results = []

println("Benchmarking across different problem sizes...")
println("=" ^ 70)
@printf("%-8s %-12s %-12s %-12s %-12s %-12s\n",
        "Cutoff", "Plane Waves", "Matrix Size", "Dense (s)", "LOBPCG (s)", "Speedup")
println("-" ^ 70)

for cutoff in cutoffs
    # Dense solver
    solver_dense = Solver(PSVWave(), geo, (64, 64); cutoff=cutoff)
    num_pw = solver_dense.basis.num_pw
    matrix_size = 2 * num_pw  # P-SV mode doubles the size

    t_dense = @elapsed begin
        bands_dense = compute_bands(solver_dense, kpath; bands=1:nbands, verbose=false)
    end

    # LOBPCG solver
    solver_lobpcg = Solver(PSVWave(), geo, (64, 64), LOBPCGMethod(); cutoff=cutoff)

    t_lobpcg = @elapsed begin
        bands_lobpcg = compute_bands(solver_lobpcg, kpath; bands=1:nbands, verbose=false)
    end

    # Accuracy check
    max_diff = maximum(abs.(bands_dense.frequencies - bands_lobpcg.frequencies))

    speedup = t_dense / t_lobpcg

    push!(results, (
        cutoff=cutoff,
        num_pw=num_pw,
        matrix_size=matrix_size,
        t_dense=t_dense,
        t_lobpcg=t_lobpcg,
        speedup=speedup,
        max_diff=max_diff,
        bands_dense=bands_dense,
        bands_lobpcg=bands_lobpcg
    ))

    @printf("%-8d %-12d %-12d %-12.2f %-12.2f %-12.2fx\n",
            cutoff, num_pw, matrix_size, t_dense, t_lobpcg, speedup)
end

println("=" ^ 70)

# ============================================================================
# Summary statistics
# ============================================================================
println("\n=== Summary ===")
avg_speedup = sum(r.speedup for r in results) / length(results)
max_error = maximum(r.max_diff for r in results)

println("Average speedup: $(round(avg_speedup, digits=2))x")
println("Maximum frequency error: $(round(max_error, digits=1)) rad/s")
println()
println("Recommendation:")
println("  - Small problems (cutoff < 10): Dense is sufficient")
println("  - Large problems (cutoff ≥ 15): LOBPCG recommended")

# ============================================================================
# Create plots
# ============================================================================

# Plot 1: Performance scaling
p1 = plot(
    xlabel="Matrix Size",
    ylabel="Time (s)",
    title="Solver Performance Scaling",
    legend=:topleft,
    yscale=:log10,
    xscale=:log10,
    size=(500, 400)
)

matrix_sizes = [r.matrix_size for r in results]
times_dense = [r.t_dense for r in results]
times_lobpcg = [r.t_lobpcg for r in results]

plot!(p1, matrix_sizes, times_dense, marker=:circle, label="Dense", linewidth=2)
plot!(p1, matrix_sizes, times_lobpcg, marker=:square, label="LOBPCG", linewidth=2)

# Plot 2: Speedup vs problem size
p2 = plot(
    xlabel="Matrix Size",
    ylabel="Speedup (Dense/LOBPCG)",
    title="LOBPCG Speedup vs Problem Size",
    legend=false,
    size=(500, 400)
)

speedups = [r.speedup for r in results]
plot!(p2, matrix_sizes, speedups, marker=:diamond, linewidth=2, color=:green)
hline!(p2, [1.0], linestyle=:dash, color=:gray, alpha=0.5)

# Plot 3: Band structure comparison (largest problem)
last_result = results[end]
p3 = plot(
    xlabel="Wave Vector",
    ylabel="ω (rad/s)",
    title="Band Structure (cutoff=$(last_result.cutoff))",
    legend=:topright,
    size=(600, 400)
)

dists = last_result.bands_dense.distances
for b in 1:size(last_result.bands_dense.frequencies, 2)
    if b == 1
        plot!(p3, dists, last_result.bands_dense.frequencies[:, b],
              color=:blue, linewidth=2, label="Dense")
        plot!(p3, dists, last_result.bands_lobpcg.frequencies[:, b],
              color=:red, linestyle=:dash, linewidth=1.5, label="LOBPCG")
    else
        plot!(p3, dists, last_result.bands_dense.frequencies[:, b],
              color=:blue, linewidth=2, label="")
        plot!(p3, dists, last_result.bands_lobpcg.frequencies[:, b],
              color=:red, linestyle=:dash, linewidth=1.5, label="")
    end
end

# Plot 4: Frequency error distribution
p4 = plot(
    xlabel="Cutoff",
    ylabel="Max Error (rad/s)",
    title="LOBPCG Accuracy vs Problem Size",
    legend=false,
    size=(500, 400)
)

cutoff_vals = [r.cutoff for r in results]
errors = [r.max_diff for r in results]
bar!(p4, cutoff_vals, errors, color=:orange, alpha=0.7)

# Combined plot
p_combined = plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))

savefig(p_combined, joinpath(@__DIR__, "208_solver_comparison.png"))
println("\nSaved: 208_solver_comparison.png")

# Also save individual performance plot
savefig(p1, joinpath(@__DIR__, "208_solver_performance.png"))
println("Saved: 208_solver_performance.png")
