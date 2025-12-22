# Last-Modified: 2025-12-18T21:00:00+09:00
# Warm start benchmark for LOBPCG solver
# Compares: random init vs low-order Fourier modes + warm start

using PhoXonic
using LinearAlgebra
using IterativeSolvers
using Printf
using Plots
using Statistics

println("=== LOBPCG Warm Start Benchmark ===")
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

steel = IsotropicElastic(; ρ=ρ_steel, λ=λ_steel, μ=μ_steel)
epoxy = IsotropicElastic(; ρ=ρ_epoxy, λ=λ_epoxy, μ=μ_epoxy)

a = 1.0
lat = hexagonal_lattice(a)
f_target = 0.4
r = sqrt(f_target * a^2 * sqrt(3) / (2π))

geo = Geometry(lat, epoxy, [(Circle([0.0, 0.0], r), steel)])

# ============================================================================
# Solver setup
# ============================================================================
cutoff = 20
solver = Solver(PSVWave(), geo, (64, 64); cutoff=cutoff)
num_pw = solver.basis.num_pw
dim = 2 * num_pw  # P-SV doubles the matrix size

println("Cutoff: $cutoff")
println("Plane waves: $num_pw")
println("Matrix size: $dim × $dim")

# k-path
kpath = simple_kpath_hexagonal(; a=a, npoints=30)
k_points = kpath.points
n_kpts = length(k_points)
println("k-points: $n_kpts")

nbands = 15
tol = 1e-4
maxiter = 300

println("Bands: $nbands")
println("Tolerance: $tol")
println()

# ============================================================================
# Helper: Create low-order Fourier mode initial vectors
# ============================================================================
function create_loworder_init(solver, nev)
    G_vecs = solver.basis.G
    num_pw = solver.basis.num_pw
    dim = 2 * num_pw  # P-SV mode

    # Sort G vectors by magnitude
    G_norms = [norm(G) for G in G_vecs]
    sorted_idx = sortperm(G_norms)

    # Create initial vectors: use lowest |G| modes
    # For P-SV, we have 2 components per G vector
    X0 = zeros(ComplexF64, dim, nev)

    for j in 1:nev
        # Alternate between x and y components
        g_idx = div(j - 1, 2) + 1
        comp = mod(j - 1, 2)  # 0 for x, 1 for y

        if g_idx <= num_pw
            pw_idx = sorted_idx[g_idx]
            X0[pw_idx + comp * num_pw, j] = 1.0
        else
            # Fallback to random if not enough G vectors
            X0[:, j] = randn(ComplexF64, dim)
        end
    end

    # Orthonormalize
    X0, _ = qr(X0)
    return Matrix(X0)
end

# ============================================================================
# Helper: Create random initial vectors
# ============================================================================
function create_random_init(dim, nev)
    X0 = randn(ComplexF64, dim, nev)
    X0, _ = qr(X0)
    return Matrix(X0)
end

# ============================================================================
# Check matrix scaling
# ============================================================================
println("--- Matrix Scaling Analysis ---")
k_test = k_points[1]
LHS_test, RHS_test = PhoXonic.build_matrices(solver, k_test)

# Characteristic scales
scale_A = maximum(abs.(LHS_test))
scale_B = maximum(abs.(RHS_test))
println("Max |A|: $(scale_A)")
println("Max |B|: $(scale_B)")
println("Ratio A/B: $(scale_A / scale_B)")

# Scaling factors to make matrices O(1)
# A_scaled = A / scale_A, B_scaled = B / scale_B
# Eigenvalue scaling: λ_original = λ_scaled * (scale_A / scale_B)
scale_factor = sqrt(scale_A / scale_B)  # For ω scaling
println("Scale factor for ω: $(scale_factor)")
println()

# ============================================================================
# Benchmark: Random initialization (no scaling)
# ============================================================================
println("--- Method 1: Random Init (no scaling) ---")

freqs_random = zeros(n_kpts, nbands)
iters_random = zeros(Int, n_kpts)
residuals_random = zeros(n_kpts)

t_random = @elapsed begin
    for (i, k) in enumerate(k_points)
        LHS, RHS = PhoXonic.build_matrices(solver, k)
        A = Hermitian(LHS)
        B = Hermitian(RHS)

        # Random initial vectors
        X0 = create_random_init(dim, nbands)

        results = IterativeSolvers.lobpcg(A, B, false, X0; tol=tol, maxiter=maxiter)

        λ_vals = results.λ
        ω² = real.(λ_vals)
        ω² = max.(ω², 0.0)
        frequencies = sqrt.(ω²)

        freqs_random[i, :] = frequencies[1:nbands]
        iters_random[i] = results.iterations
        residuals_random[i] = maximum(results.residual_norms)
    end
end

println("Time: $(round(t_random, digits=2)) s")
println("Avg iterations: $(round(mean(iters_random), digits=1))")
println("Max residual: $(round(maximum(residuals_random), sigdigits=3))")

# ============================================================================
# Benchmark: Scaled matrices + Warm start + Diagonal Preconditioner
# ============================================================================
println("\n--- Method 2: Scaled + Warm Start + Preconditioner ---")

freqs_warmstart = zeros(n_kpts, nbands)
iters_warmstart = zeros(Int, n_kpts)
residuals_warmstart = zeros(n_kpts)
prev_eigvecs = nothing

t_warmstart = @elapsed begin
    for (i, k) in enumerate(k_points)
        LHS, RHS = PhoXonic.build_matrices(solver, k)

        if i == 1
            # First k-point: use Dense for accurate initial solution
            λ_all, V_all = eigen(Hermitian(LHS), Hermitian(RHS))
            ω² = real.(λ_all)
            ω² = max.(ω², 0.0)
            frequencies = sqrt.(ω²)

            freqs_warmstart[i, :] = frequencies[1:nbands]
            iters_warmstart[i] = 0  # Dense = 0 iterations
            residuals_warmstart[i] = 0.0

            # Store eigenvectors for warm start
            global prev_eigvecs = V_all[:, 1:nbands]
        else
            # Subsequent k-points: use LOBPCG with warm start
            # Scale A only (to preserve B positive definiteness)
            A_scaled = Hermitian(LHS / scale_A)
            B_scaled = Hermitian(RHS)

            # Diagonal preconditioner on scaled matrix
            d = diag(A_scaled)
            d_safe = [abs(x) > 1e-10 ? x : 1e-10 for x in d]
            P = Diagonal(1.0 ./ d_safe)

            # Use previous eigenvectors as initial guess
            X0 = prev_eigvecs
            X0, _ = qr(X0)
            X0 = Matrix(X0)

            results = IterativeSolvers.lobpcg(
                A_scaled, B_scaled, false, X0; P=P, tol=tol, maxiter=maxiter
            )

            # Rescale eigenvalues back to original
            λ_scaled = results.λ
            λ_original = λ_scaled * scale_A
            ω² = real.(λ_original)
            ω² = max.(ω², 0.0)
            frequencies = sqrt.(ω²)

            freqs_warmstart[i, :] = frequencies[1:nbands]
            iters_warmstart[i] = results.iterations
            residuals_warmstart[i] = maximum(results.residual_norms)

            # Store eigenvectors for next k-point
            global prev_eigvecs = results.X
        end
    end
end

println("Time: $(round(t_warmstart, digits=2)) s")
println("Avg iterations: $(round(mean(iters_warmstart), digits=1))")
println("Max residual: $(round(maximum(residuals_warmstart), sigdigits=3))")

# ============================================================================
# Comparison with Dense solver
# ============================================================================
println("\n--- Method 3: Dense (Reference) ---")

freqs_dense = zeros(n_kpts, nbands)

t_dense = @elapsed begin
    for (i, k) in enumerate(k_points)
        LHS, RHS = PhoXonic.build_matrices(solver, k)

        λ_all, V_all = eigen(Hermitian(LHS), Hermitian(RHS))

        ω² = real.(λ_all)
        ω² = max.(ω², 0.0)
        frequencies = sqrt.(ω²)

        freqs_dense[i, :] = frequencies[1:nbands]
    end
end

println("Time: $(round(t_dense, digits=2)) s")

# ============================================================================
# Results Summary
# ============================================================================
println("\n" * "=" ^ 60)
println("SUMMARY")
println("=" ^ 60)

# Accuracy comparison (vs Dense)
error_random = maximum(abs.(freqs_random - freqs_dense))
error_warmstart = maximum(abs.(freqs_warmstart - freqs_dense))

@printf("%-25s %12s %12s %12s\n", "Method", "Time (s)", "Speedup", "Max Error")
println("-" ^ 60)
@printf("%-25s %12.2f %12s %12s\n", "Dense (reference)", t_dense, "1.0x", "-")
@printf(
    "%-25s %12.2f %12.2fx %12.1f\n",
    "LOBPCG (random)",
    t_random,
    t_dense/t_random,
    error_random
)
@printf(
    "%-25s %12.2f %12.2fx %12.1f\n",
    "LOBPCG (warm start)",
    t_warmstart,
    t_dense/t_warmstart,
    error_warmstart
)
println("=" ^ 60)

# ============================================================================
# Plot: Iteration count comparison
# ============================================================================
p1 = plot(;
    xlabel="k-point index",
    ylabel="Iterations",
    title="LOBPCG Iterations per k-point",
    legend=:topright,
    size=(600, 400),
)

plot!(p1, 1:n_kpts, iters_random; label="Random init", marker=:circle, markersize=3)
plot!(p1, 1:n_kpts, iters_warmstart; label="Warm start", marker=:square, markersize=3)

# ============================================================================
# Plot: Band structure comparison
# ============================================================================
dists = kpath.distances

p2 = plot(;
    xlabel="Wave Vector",
    ylabel="ω (rad/s)",
    title="Band Structure Comparison",
    legend=:topright,
    size=(600, 400),
)

for b in 1:nbands
    if b == 1
        plot!(p2, dists, freqs_dense[:, b]; color=:black, linewidth=2, label="Dense")
        plot!(
            p2,
            dists,
            freqs_warmstart[:, b];
            color=:red,
            linestyle=:dash,
            label="Warm start",
        )
    else
        plot!(p2, dists, freqs_dense[:, b]; color=:black, linewidth=2, label="")
        plot!(p2, dists, freqs_warmstart[:, b]; color=:red, linestyle=:dash, label="")
    end
end

# ============================================================================
# Combined plot
# ============================================================================
p_combined = plot(p1, p2; layout=(1, 2), size=(1100, 400))

savefig(p_combined, joinpath(@__DIR__, "209_warmstart_benchmark.png"))
println("\nSaved: 209_warmstart_benchmark.png")

# ============================================================================
# Detailed iteration statistics
# ============================================================================
println("\nIteration Statistics:")
println("  Random init:")
println("    First k-point: $(iters_random[1])")
println(
    "    Min: $(minimum(iters_random)), Max: $(maximum(iters_random)), Avg: $(round(mean(iters_random), digits=1))",
)
println("  Warm start:")
println("    First k-point: $(iters_warmstart[1])")
println(
    "    Min: $(minimum(iters_warmstart)), Max: $(maximum(iters_warmstart)), Avg: $(round(mean(iters_warmstart), digits=1))",
)
