# Last-Modified: 2025-12-25T15:00:00+09:00
# Detailed solver comparison: Dense vs LOBPCG
# Compares performance and accuracy across different problem sizes
#
# Key findings:
# - LOBPCG becomes faster than Dense for cutoff ≥ 12
# - warm_start=true (default) is essential for good performance
# - Typical accuracy: 10-100 rad/s error (0.1-1% relative)

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

println("=== Detailed Solver Comparison ===")
println("Dense vs LOBPCG performance analysis")
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

# k-path with moderate resolution
kpath = simple_kpath_hexagonal(; a=a, npoints=20)

# ============================================================================
# Benchmark across different cutoff values
# ============================================================================
cutoffs = [8, 10, 12, 15]
nbands = 10

results = []

println("Benchmarking across different problem sizes...")
println("=" ^ 75)
@printf(
    "%-8s %-10s %-10s %-10s %-10s %-10s %-12s\n",
    "Cutoff",
    "N(pw)",
    "Matrix",
    "Dense(s)",
    "LOBPCG(s)",
    "Speedup",
    "Error(rad/s)"
)
println("-" ^ 75)

for cutoff in cutoffs
    # Dense solver
    solver_dense = Solver(PSVWave(), geo, (64, 64); cutoff=cutoff)
    num_pw = solver_dense.basis.num_pw
    matrix_size = 2 * num_pw  # P-SV mode doubles the size

    t_dense = @elapsed begin
        bands_dense = compute_bands(solver_dense, kpath; bands=1:nbands, verbose=false)
    end

    # LOBPCG solver with default settings
    solver_lobpcg = Solver(PSVWave(), geo, (64, 64), LOBPCGMethod(); cutoff=cutoff)

    t_lobpcg = @elapsed begin
        bands_lobpcg = compute_bands(solver_lobpcg, kpath; bands=1:nbands, verbose=false)
    end

    # Accuracy check
    max_diff = maximum(abs.(bands_dense.frequencies - bands_lobpcg.frequencies))

    speedup = t_dense / t_lobpcg

    push!(
        results,
        (
            cutoff=cutoff,
            num_pw=num_pw,
            matrix_size=matrix_size,
            t_dense=t_dense,
            t_lobpcg=t_lobpcg,
            speedup=speedup,
            max_diff=max_diff,
            bands_dense=bands_dense,
            bands_lobpcg=bands_lobpcg,
        ),
    )

    speedup_str = speedup >= 1 ? @sprintf("%.2fx", speedup) : @sprintf("%.2fx*", speedup)
    @printf(
        "%-8d %-10d %-10d %-10.1f %-10.1f %-10s %-12.1f\n",
        cutoff,
        num_pw,
        matrix_size,
        t_dense,
        t_lobpcg,
        speedup_str,
        max_diff
    )
end

println("=" ^ 75)
println("* speedup < 1 means Dense is faster")

# ============================================================================
# Summary statistics
# ============================================================================
println("\n=== Summary ===")
max_error = maximum(r.max_diff for r in results)
best_speedup_idx = argmax([r.speedup for r in results])
best_result = results[best_speedup_idx]

println("Maximum frequency error: $(round(max_error, digits=1)) rad/s")
println(
    "Best speedup: $(round(best_result.speedup, digits=2))x at cutoff=$(best_result.cutoff)",
)
println()
println("Recommendation:")
println("  - cutoff < 10: Dense is sufficient")
println("  - cutoff ≥ 12: LOBPCG recommended (speedup increases with problem size)")

# ============================================================================
# Create plots
# ============================================================================

matrix_sizes = [r.matrix_size for r in results]
times_dense = [r.t_dense for r in results]
times_lobpcg = [r.t_lobpcg for r in results]
speedups = [r.speedup for r in results]
cutoff_vals = [r.cutoff for r in results]
errors = [r.max_diff for r in results]

# Plot 1: Time comparison (log scale)
p1 = plot(;
    xlabel="Matrix Size",
    ylabel="Time (s)",
    title="Solver Performance",
    legend=:topleft,
    yscale=:log10,
    size=(600, 450),
)
plot!(
    p1, matrix_sizes, times_dense; marker=:circle, label="Dense", linewidth=2, color=:blue
)
plot!(
    p1,
    matrix_sizes,
    times_lobpcg;
    marker=:square,
    label="LOBPCG",
    linewidth=2,
    color=:red,
)

# Plot 2: Speedup
p2 = plot(;
    xlabel="Cutoff",
    ylabel="Speedup (Dense/LOBPCG)",
    title="LOBPCG Speedup",
    legend=false,
    size=(600, 450),
)
bar!(p2, cutoff_vals, speedups; color=:green, alpha=0.7)
hline!(p2, [1.0]; linestyle=:dash, color=:gray, linewidth=2, label="")
annotate!(p2, [(cutoff_vals[end], 1.1, text("Dense = LOBPCG", 10, :gray))])

# Plot 3: Accuracy
p3 = plot(;
    xlabel="Cutoff",
    ylabel="Max Error (rad/s)",
    title="Frequency Error",
    legend=false,
    size=(600, 450),
)
bar!(p3, cutoff_vals, errors; color=:orange, alpha=0.7)

# Plot 4: Band structure comparison (largest cutoff)
last_result = results[end]
dists = last_result.bands_dense.distances
labels_data = last_result.bands_dense.labels
label_positions = [dists[i] for (i, _) in labels_data]
label_names = [l for (_, l) in labels_data]

# Convert to krad/s
freq_dense_k = last_result.bands_dense.frequencies / 1000
freq_lobpcg_k = last_result.bands_lobpcg.frequencies / 1000

p4 = plot(;
    xlabel="Wave Vector",
    ylabel="ω (10³ rad/s)",
    title="Band Structure (cutoff=$(last_result.cutoff))",
    legend=:topright,
    size=(600, 450),
)

# Plot Dense bands
for b in 1:size(freq_dense_k, 2)
    label = b == 1 ? "Dense" : ""
    plot!(p4, dists, freq_dense_k[:, b]; color=:blue, linewidth=2, label=label)
end

# Plot LOBPCG bands (dashed)
for b in 1:size(freq_lobpcg_k, 2)
    label = b == 1 ? "LOBPCG" : ""
    plot!(
        p4,
        dists,
        freq_lobpcg_k[:, b];
        color=:red,
        linewidth=1,
        linestyle=:dash,
        label=label,
    )
end

vline!(p4, label_positions; color=:gray, linestyle=:dot, alpha=0.5, label="")
xticks!(p4, label_positions, label_names)

# Combined plot: 2x2 layout
p_combined = plot(
    p1,
    p2,
    p3,
    p4;
    layout=(2, 2),
    size=(1200, 900),
    plot_title="Solver Comparison: Steel/Epoxy PSV Wave",
)

savefig(p_combined, joinpath(@__DIR__, "208_solver_comparison.png"))
println("\nSaved: 208_solver_comparison.png")
