# Last-Modified: 2025-12-15T21:30:00+09:00
# PhoXonic.jl Extension for ReducedShiftedKrylov.jl
# This extension provides RSKGF and MatrixFreeGF methods for Green's function computation.

module PhoXonicReducedShiftedKrylovExt

using PhoXonic
using ReducedShiftedKrylov
using LinearAlgebra: dot, norm

# Import types from PhoXonic
import PhoXonic: Solver, Dim1, Dim2, Dim3, WaveType
import PhoXonic: GFMethod, RSKGF, MatrixFreeGF, DirectGF
import PhoXonic: RHSInvMethod, ApproximateRHSInv, CGRHSInv
import PhoXonic: EffectiveHamiltonian, MatrixFreeEffectiveHamiltonian, NegatedOperator
import PhoXonic: MatrixFreeOperator, FFTContext, MatrixFreeWorkspace
import PhoXonic: build_matrices, _to_kvec, _to_kvec_f64

# Import functions to extend
import PhoXonic: compute_greens_function, compute_dos, compute_ldos

# ============================================================================
# compute_greens_function implementations
# ============================================================================

function PhoXonic.compute_greens_function(solver::Solver, k, ω_values::AbstractVector{<:Real},
                                          source::AbstractVector, method::RSKGF; η::Real=1e-3)
    # Build effective Hamiltonian: H = RHS⁻¹ * LHS
    H = EffectiveHamiltonian(solver, k)

    # Negate for RSCG interface: A = -H
    # RSCG solves (σI + A)x = b, so with A = -H:
    # (σI - H)x = b  which is what we want for G(z) = (z - H)⁻¹
    A = NegatedOperator(H)

    # Shifts: σ_j = ω_j² + iη
    shifts = ComplexF64[ω^2 + im*η for ω in ω_values]

    # Solve using ReducedShiftedKrylov.jl
    x_solutions, _ = ReducedShiftedKrylov.rscg(A, ComplexF64.(source), shifts;
                                                atol=method.atol, rtol=method.rtol,
                                                itmax=method.itmax, verbose=method.verbose)

    return x_solutions
end

function PhoXonic.compute_greens_function(solver::Solver{Dim2,W}, k, ω_values::AbstractVector{<:Real},
                                          source::AbstractVector, method::MatrixFreeGF;
                                          η::Real=1e-3) where W<:WaveType
    k_vec = _to_kvec_f64(Dim2, k)

    # Create FFT context and workspace
    ctx = FFTContext(solver)
    workspace = MatrixFreeWorkspace(ctx)

    # Create matrix-free operator
    op = MatrixFreeOperator(solver, Vector{Float64}(k_vec), ctx, workspace)

    # Create matrix-free effective Hamiltonian
    H = MatrixFreeEffectiveHamiltonian(op, method.rhs_inv_method)

    # Negate for RSCG interface: A = -H
    A = NegatedOperator(H)

    # Shifts: σ_j = ω_j² + iη
    shifts = ComplexF64[ω^2 + im*η for ω in ω_values]

    # Solve using ReducedShiftedKrylov.jl
    x_solutions, _ = ReducedShiftedKrylov.rscg(A, ComplexF64.(source), shifts;
                                                atol=method.atol, rtol=method.rtol,
                                                itmax=method.itmax, verbose=method.verbose)

    return x_solutions
end

function PhoXonic.compute_greens_function(solver::Solver{Dim3,W}, k, ω_values::AbstractVector{<:Real},
                                          source::AbstractVector, method::MatrixFreeGF;
                                          η::Real=1e-3) where W<:WaveType
    k_vec = _to_kvec_f64(Dim3, k)

    # Create FFT context and workspace
    ctx = FFTContext(solver)
    workspace = MatrixFreeWorkspace(ctx)

    # Create matrix-free operator
    op = MatrixFreeOperator(solver, Vector{Float64}(k_vec), ctx, workspace)

    # Create matrix-free effective Hamiltonian
    H = MatrixFreeEffectiveHamiltonian(op, method.rhs_inv_method)

    # Negate for RSCG interface: A = -H
    A = NegatedOperator(H)

    # Shifts: σ_j = ω_j² + iη
    shifts = ComplexF64[ω^2 + im*η for ω in ω_values]

    # Solve using ReducedShiftedKrylov.jl
    x_solutions, _ = ReducedShiftedKrylov.rscg(A, ComplexF64.(source), shifts;
                                                atol=method.atol, rtol=method.rtol,
                                                itmax=method.itmax, verbose=method.verbose)

    return x_solutions
end

# ============================================================================
# compute_dos implementations (stochastic trace for RSKGF/MatrixFreeGF)
# ============================================================================

function PhoXonic.compute_dos(solver::Solver{Dim2}, ω_values::AbstractVector{<:Real},
                              k_points::AbstractVector, method::RSKGF;
                              η::Real=1e-3, n_random::Int=10)
    N = solver.basis.num_pw
    dos = zeros(length(ω_values))

    for k in k_points
        for _ in 1:n_random
            # Random source vector
            source = randn(ComplexF64, N)
            source ./= norm(source)

            G_values = compute_greens_function(solver, k, ω_values, source, method; η=η)

            for (iω, G) in enumerate(G_values)
                dos[iω] += -imag(dot(source, G)) / π
            end
        end
    end

    dos ./= (length(k_points) * n_random)
    return dos
end

function PhoXonic.compute_dos(solver::Solver{Dim2}, ω_values::AbstractVector{<:Real},
                              k_points::AbstractVector, method::MatrixFreeGF;
                              η::Real=1e-3, n_random::Int=10)
    N = solver.basis.num_pw
    dos = zeros(length(ω_values))

    for k in k_points
        for _ in 1:n_random
            # Random source vector
            source = randn(ComplexF64, N)
            source ./= norm(source)

            G_values = compute_greens_function(solver, k, ω_values, source, method; η=η)

            for (iω, G) in enumerate(G_values)
                dos[iω] += -imag(dot(source, G)) / π
            end
        end
    end

    dos ./= (length(k_points) * n_random)
    return dos
end

# ============================================================================
# compute_ldos implementations
# ============================================================================

function PhoXonic.compute_ldos(solver::Solver{Dim2}, position::AbstractVector{<:Real},
                               ω_values::AbstractVector{<:Real},
                               k_points::AbstractVector,
                               method::RSKGF;
                               η::Real=1e-3)
    N = solver.basis.num_pw
    ldos = zeros(length(ω_values))

    # Source at given position
    r = position
    source = ComplexF64[exp(-im * dot(G, r)) for G in solver.basis.G]
    source ./= N

    for k in k_points
        k_vec = _to_kvec_f64(Dim2, k)

        # Compute Green's function using RSKGF method
        G_values = compute_greens_function(solver, k_vec, ω_values, source, method; η=η)

        for (iω, G) in enumerate(G_values)
            G_rr = dot(conj.(source), G)
            ldos[iω] += -imag(G_rr) / π
        end
    end

    ldos ./= length(k_points)
    return ldos
end

function PhoXonic.compute_ldos(solver::Solver{Dim2,W}, position::AbstractVector{<:Real},
                               ω_values::AbstractVector{<:Real},
                               k_points::AbstractVector,
                               method::MatrixFreeGF;
                               η::Real=1e-3) where W<:WaveType
    N = solver.basis.num_pw
    ldos = zeros(length(ω_values))

    # Source at given position
    r = position
    source = ComplexF64[exp(-im * dot(G, r)) for G in solver.basis.G]
    source ./= N

    # Create FFT context and workspace once, reuse for all k-points
    ctx = FFTContext(solver)
    workspace = MatrixFreeWorkspace(ctx)

    for k in k_points
        k_vec = _to_kvec_f64(Dim2, k)

        # Create matrix-free operator (reuse ctx and workspace)
        op = MatrixFreeOperator(solver, Vector{Float64}(k_vec), ctx, workspace)
        H = MatrixFreeEffectiveHamiltonian(op, method.rhs_inv_method)
        A = NegatedOperator(H)

        shifts = ComplexF64[ω^2 + im*η for ω in ω_values]
        G_values, _ = ReducedShiftedKrylov.rscg(A, ComplexF64.(source), shifts;
                                                 atol=method.atol, rtol=method.rtol,
                                                 itmax=method.itmax, verbose=method.verbose)

        for (iω, G) in enumerate(G_values)
            G_rr = dot(conj.(source), G)
            ldos[iω] += -imag(G_rr) / π
        end
    end

    ldos ./= length(k_points)
    return ldos
end

# 3D support for MatrixFreeGF
function PhoXonic.compute_ldos(solver::Solver{Dim3,W}, position::AbstractVector{<:Real},
                               ω_values::AbstractVector{<:Real},
                               k_points::AbstractVector,
                               method::MatrixFreeGF;
                               η::Real=1e-3) where W<:WaveType
    N = solver.basis.num_pw
    nc = PhoXonic.ncomponents(solver.wave)  # 3 for 3D waves
    ldos = zeros(length(ω_values))

    # Source at given position (for each component)
    r = position
    source_single = ComplexF64[exp(-im * dot(G, r)) for G in solver.basis.G]
    source_single ./= N

    # Full source vector (all components)
    source = zeros(ComplexF64, N * nc)
    for c in 1:nc
        source[(c-1)*N+1:c*N] .= source_single
    end

    # Create FFT context and workspace once, reuse for all k-points
    ctx = FFTContext(solver)
    workspace = MatrixFreeWorkspace(ctx)

    for k in k_points
        k_vec = _to_kvec_f64(Dim3, k)

        # Create matrix-free operator (reuse ctx and workspace)
        op = MatrixFreeOperator(solver, Vector{Float64}(k_vec), ctx, workspace)
        H = MatrixFreeEffectiveHamiltonian(op, method.rhs_inv_method)
        A = NegatedOperator(H)

        shifts = ComplexF64[ω^2 + im*η for ω in ω_values]
        G_values, _ = ReducedShiftedKrylov.rscg(A, ComplexF64.(source), shifts;
                                                 atol=method.atol, rtol=method.rtol,
                                                 itmax=method.itmax, verbose=method.verbose)

        for (iω, G) in enumerate(G_values)
            G_rr = dot(conj.(source), G)
            ldos[iω] += -imag(G_rr) / π
        end
    end

    ldos ./= length(k_points)
    return ldos
end

end # module
