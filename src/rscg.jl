# Last-Modified: 2025-12-15T22:00:00+09:00

#=
Reduced Shifted Conjugate Gradient (RSCG) method for Green's function computation.

The RSCG method efficiently solves multiple shifted linear systems:
    (z_j - A) x_j = b  for j = 1, 2, ..., m

where A is Hermitian and z_j are complex shifts (ω² + iη in our case).

For photonic/phononic crystals, we solve:
    G(ω) = (ω² + iη - H)⁻¹
where H = RHS⁻¹ * LHS is the effective Hamiltonian.

References:
- Y. Nagai et al., J. Phys. Soc. Jpn. 86, 014708 (2017)
- Y. Nagai et al., Springer (2017) "Reduced-Shifted Conjugate-Gradient Method for a Green's Function"
=#

# ============================================================================
# K-point conversion helpers (multiple dispatch instead of isa checks)
# ============================================================================

"""
    _to_kvec(k) -> Vector

Convert k-point to Vector. Uses multiple dispatch instead of runtime type checks.
"""
_to_kvec(k::AbstractVector) = k
_to_kvec(k) = [k[1], k[2]]  # For tuples, named tuples, BrillouinZoneCoordinate, etc.

"""
    _to_kvec_f64(::Type{D}, k) -> Vector{Float64}

Convert k-point to Vector{Float64} for the given dimension.
"""
_to_kvec_f64(::Type{Dim2}, k::AbstractVector) = Vector{Float64}(k)
_to_kvec_f64(::Type{Dim2}, k) = Float64[k[1], k[2]]
_to_kvec_f64(::Type{Dim3}, k::AbstractVector) = Vector{Float64}(k)
_to_kvec_f64(::Type{Dim3}, k) = Float64[k[1], k[2], k[3]]

"""
    _apply_operator(A, p) -> Vector

Apply linear operator A to vector p. Uses multiple dispatch for functions vs matrices.
"""
_apply_operator(A::Function, p) = A(p)
_apply_operator(A, p) = A * p  # Fallback for matrices/AbstractMatrix

# ============================================================================
# RHS⁻¹ Method Types (for multiple dispatch)
# ============================================================================

"""
    RHSInvMethod

Abstract type for RHS⁻¹ application methods in matrix-free effective Hamiltonian.

Concrete subtypes:
- [`ApproximateRHSInv`](@ref): Fast element-wise 1/ε in real space (approximate)
- [`CGRHSInv`](@ref): Iterative CG to solve RHS * y = x (exact)

See also: [`MatrixFreeEffectiveHamiltonian`](@ref), [`MatrixFreeGF`](@ref)
"""
abstract type RHSInvMethod end

"""
    ApproximateRHSInv <: RHSInvMethod

Approximate method for RHS⁻¹ application: element-wise 1/ε (or 1/ρ) in real space.

Fast but inaccurate for high-contrast media (e.g., Si/Air with ε ratio ≈ 12:1).

# Example
```julia
method = MatrixFreeGF(rhs_inv_method=ApproximateRHSInv())
```
"""
struct ApproximateRHSInv <: RHSInvMethod end

"""
    CGRHSInv <: RHSInvMethod

CG-based method for RHS⁻¹ application: solve RHS * y = x iteratively.

Accurate but slower than approximate method. Recommended for high-contrast media.

# Fields
- `atol::Float64`: Absolute tolerance (default: 1e-10)
- `rtol::Float64`: Relative tolerance (default: 1e-10)
- `maxiter::Int`: Maximum iterations (default: 100)

# Example
```julia
# Default parameters
method = MatrixFreeGF(rhs_inv_method=CGRHSInv())

# Custom parameters for high-contrast media
method = MatrixFreeGF(rhs_inv_method=CGRHSInv(atol=1e-12, rtol=1e-12, maxiter=200))
```
"""
struct CGRHSInv <: RHSInvMethod
    atol::Float64
    rtol::Float64
    maxiter::Int
end
function CGRHSInv(; atol::Float64=1e-10, rtol::Float64=1e-10, maxiter::Int=100)
    CGRHSInv(atol, rtol, maxiter)
end

# Conversion from Symbol (deprecated, for backward compatibility)
function _convert_rhs_inv_method(method::Symbol)
    Base.depwarn(
        "Symbol-based rhs_inv_method is deprecated. Use ApproximateRHSInv() or CGRHSInv() instead.",
        :_convert_rhs_inv_method,
    )
    if method == :approximate
        return ApproximateRHSInv()
    elseif method == :cg
        return CGRHSInv()
    else
        throw(ArgumentError("Unknown rhs_inv_method: $method. Use :approximate or :cg"))
    end
end
_convert_rhs_inv_method(method::RHSInvMethod) = method

# Error fallback for unsupported types
function _convert_rhs_inv_method(method)
    throw(
        ArgumentError(
            "Unsupported rhs_inv_method type: $(typeof(method)). Use ApproximateRHSInv() or CGRHSInv().",
        ),
    )
end

"""
    compute_greens_function(solver::Solver, k, ω_values, source; η=1e-3)

Compute Green's function G(ω) = (ω² + iη - H)⁻¹ * source for multiple frequencies.

Uses direct solve (dense method). For large systems, use matrix-free iterative methods.

# Arguments
- `solver`: The PhoXonic solver
- `k`: Wave vector
- `ω_values`: Frequencies at which to compute G
- `source`: Source vector in plane wave basis
- `η`: Broadening parameter (small positive imaginary part)

# Returns
- Vector of Green's function values G(ω) * source for each ω
"""
function compute_greens_function(
    solver::Solver,
    k,
    ω_values::AbstractVector{<:Real},
    source::AbstractVector;
    η::Real=1e-3,
)
    N = solver.basis.num_pw
    T = ComplexF64

    # Build matrices: LHS * ψ = ω² * RHS * ψ
    # Effective Hamiltonian: H = RHS⁻¹ * LHS
    # Green's function: G(ω) = (ω² + iη - H)⁻¹
    # We solve: (ω² + iη) RHS - LHS) x = RHS * source

    LHS, RHS = build_matrices(solver, k)
    RHS_source = RHS * source

    results = Vector{Vector{T}}(undef, length(ω_values))

    for (i, ω) in enumerate(ω_values)
        z = ω^2 + im*η
        # (z * RHS - LHS) * x = RHS * source
        A = z * RHS - LHS
        results[i] = A \ RHS_source
    end

    return results
end

# Convenience method with MatrixFreeOperator
function compute_greens_function(
    op::MatrixFreeOperator,
    ω_values::AbstractVector{<:Real},
    source::AbstractVector;
    η::Real=1e-3,
)
    compute_greens_function(op.solver, op.k, ω_values, source; η=η)
end

# ============================================================================
# DOS and LDOS computation
# ============================================================================

"""
    compute_dos(solver::Solver, ω_values, k_points; η=1e-3)

Compute density of states (DOS) by summing over k-points.

DOS(ω) = -1/π Im[Tr G(ω)] summed over k-points

# Arguments
- `solver`: The PhoXonic solver
- `ω_values`: Frequencies at which to compute DOS
- `k_points`: List of k-points to sum over (should cover Brillouin zone)
- `η`: Broadening parameter

# Returns
- Vector of DOS values at each frequency
"""
function compute_dos(
    solver::Solver{Dim2},
    ω_values::AbstractVector{<:Real},
    k_points::AbstractVector;
    η::Real=1e-3,
)
    N = solver.basis.num_pw
    dos = zeros(length(ω_values))

    # Identity source to compute Tr G
    # Tr G = Σ_n G_nn = Σ_n <e_n | G | e_n>
    # We sum diagonal elements

    for k in k_points
        k_vec = _to_kvec(k)

        # Build matrices once per k-point
        LHS, RHS = build_matrices(solver, k_vec)

        for (iω, ω) in enumerate(ω_values)
            z = ω^2 + im*η
            A = z * RHS - LHS

            # Compute Tr G by solving for each basis vector
            # This is expensive but accurate
            # For efficiency, could use stochastic trace estimation
            A_factored = lu(A)
            trace_G = zero(ComplexF64)

            for n in 1:N
                e_n = zeros(ComplexF64, N)
                e_n[n] = 1.0
                b = RHS * e_n
                x = A_factored \ b
                trace_G += x[n]
            end

            # DOS contribution: -1/π Im[Tr G]
            dos[iω] += -imag(trace_G) / π
        end
    end

    # Normalize by number of k-points
    dos ./= length(k_points)

    return dos
end

"""
    compute_ldos(solver::Solver, position, ω_values, k_points; η=1e-3)

Compute local density of states (LDOS) at a given position.

LDOS(r, ω) = -1/π Im[G(r, r, ω)]

# Arguments
- `solver`: The PhoXonic solver
- `position`: Position in real space (Vec2 or [x, y])
- `ω_values`: Frequencies at which to compute LDOS
- `k_points`: List of k-points to sum over
- `η`: Broadening parameter

# Returns
- Vector of LDOS values at each frequency
"""
function compute_ldos(
    solver::Solver{Dim2},
    position::AbstractVector{<:Real},
    ω_values::AbstractVector{<:Real},
    k_points::AbstractVector;
    η::Real=1e-3,
)
    N = solver.basis.num_pw
    ldos = zeros(length(ω_values))

    # Source at given position in real space
    # Convert to plane wave coefficients: source[G] = e^{-iG·r}
    r = position
    source = [exp(-im * dot(G, r)) for G in solver.basis.G]
    source ./= N  # Normalization

    for k in k_points
        k_vec = _to_kvec(k)

        # Compute Green's function at this k-point
        G_values = compute_greens_function(solver, k_vec, ω_values, source; η=η)

        for (iω, G) in enumerate(G_values)
            # LDOS = -1/π Im[G(r, r)]
            # G(r, r) = Σ_G Σ_G' e^{iG·r} G_{GG'} e^{-iG'·r}
            # With our source: G(r,r) ∝ source' * G_result
            G_rr = dot(conj.(source), G)
            ldos[iω] += -imag(G_rr) / π
        end
    end

    # Normalize by number of k-points
    ldos ./= length(k_points)

    return ldos
end

# 1D versions
function compute_dos(
    solver::Solver{Dim1},
    ω_values::AbstractVector{<:Real},
    k_points::AbstractVector{<:Real};
    η::Real=1e-3,
)
    N = solver.basis.num_pw
    dos = zeros(length(ω_values))

    for k in k_points
        LHS, RHS = build_matrices(solver, k)

        for (iω, ω) in enumerate(ω_values)
            z = ω^2 + im*η
            A = z * RHS - LHS
            A_factored = lu(A)
            trace_G = zero(ComplexF64)

            for n in 1:N
                e_n = zeros(ComplexF64, N)
                e_n[n] = 1.0
                b = RHS * e_n
                x = A_factored \ b
                trace_G += x[n]
            end

            dos[iω] += -imag(trace_G) / π
        end
    end

    dos ./= length(k_points)
    return dos
end

function compute_ldos(
    solver::Solver{Dim1},
    position::Real,
    ω_values::AbstractVector{<:Real},
    k_points::AbstractVector{<:Real};
    η::Real=1e-3,
)
    N = solver.basis.num_pw
    ldos = zeros(length(ω_values))

    r = position
    source = [exp(-im * G[1] * r) for G in solver.basis.G]
    source ./= N

    for k in k_points
        G_values = compute_greens_function(solver, k, ω_values, source; η=η)

        for (iω, G) in enumerate(G_values)
            G_rr = dot(conj.(source), G)
            ldos[iω] += -imag(G_rr) / π
        end
    end

    ldos ./= length(k_points)
    return ldos
end

# ============================================================================
# Stochastic DOS computation
# ============================================================================

"""
    compute_dos_stochastic(solver::Solver, ω_values, k_points; η=1e-3, n_random=10)

Compute DOS using stochastic trace estimation.

This approximates Tr[G] by averaging over random vectors:
    Tr[G] ≈ (1/n_random) Σ_i <v_i | G | v_i>

More efficient than exact trace for large systems.

# Arguments
- `solver`: The PhoXonic solver
- `ω_values`: Frequencies at which to compute DOS
- `k_points`: List of k-points
- `η`: Broadening parameter
- `n_random`: Number of random vectors for stochastic averaging
"""
function compute_dos_stochastic(
    solver::Solver{Dim2},
    ω_values::AbstractVector{<:Real},
    k_points::AbstractVector;
    η::Real=1e-3,
    n_random::Int=10,
)
    N = solver.basis.num_pw
    dos = zeros(length(ω_values))

    for k in k_points
        k_vec = _to_kvec(k)
        LHS, RHS = build_matrices(solver, k_vec)

        # Generate random vectors (Rademacher: ±1 with equal probability)
        for _ in 1:n_random
            v = ComplexF64.(2.0 .* (rand(N) .> 0.5) .- 1.0)
            v ./= sqrt(N)

            RHS_v = RHS * v

            for (iω, ω) in enumerate(ω_values)
                z = ω^2 + im*η
                A = z * RHS - LHS
                x = A \ RHS_v
                # <v | G | v> ≈ v' * x
                trace_contribution = dot(v, x)
                dos[iω] += -imag(trace_contribution) / π
            end
        end
    end

    dos ./= (length(k_points) * n_random)
    return dos
end

# ============================================================================
# EffectiveHamiltonian: H = RHS⁻¹ * LHS
# ============================================================================

"""
    EffectiveHamiltonian{TL, TR, TF}

Effective Hamiltonian operator H = RHS⁻¹ * LHS for use with Krylov.jl.

Transforms the generalized eigenvalue problem:
    LHS * ψ = ω² * RHS * ψ
into the standard form:
    H * ψ = ω² * ψ

For Green's function computation:
    G(z) = (z·I - H)⁻¹ = (z·I - RHS⁻¹·LHS)⁻¹

# Fields
- `LHS`: Left-hand side matrix
- `RHS`: Right-hand side matrix
- `RHS_factorization`: LU factorization of RHS for efficient solves
- `n`: Problem dimension
- `tmp`: Temporary vector for intermediate computations

# Usage with Krylov.jl
Use `NegatedOperator(H)` to get `-H`, then solve `(σI + A)x = b` where `A = -H`.
"""
struct EffectiveHamiltonian{TL<:AbstractMatrix,TR<:AbstractMatrix,TF,TV<:AbstractVector}
    LHS::TL
    RHS::TR
    RHS_factorization::TF
    n::Int
    tmp::TV
end

"""
    EffectiveHamiltonian(LHS, RHS)

Create an EffectiveHamiltonian from LHS and RHS matrices.
Computes LU factorization of RHS for efficient `RHS⁻¹` application.
"""
function EffectiveHamiltonian(LHS::AbstractMatrix{T}, RHS::AbstractMatrix{T}) where {T}
    n = size(LHS, 1)
    @assert size(LHS) == size(RHS) == (n, n) "LHS and RHS must be square and same size"
    RHS_fact = lu(RHS)
    tmp = Vector{T}(undef, n)
    EffectiveHamiltonian(LHS, RHS, RHS_fact, n, tmp)
end

"""
    EffectiveHamiltonian(solver::Solver, k)

Create an EffectiveHamiltonian from a PhoXonic Solver at wave vector k.
"""
function EffectiveHamiltonian(solver::Solver, k)
    LHS, RHS = build_matrices(solver, k)
    EffectiveHamiltonian(Matrix(LHS), Matrix(RHS))
end

# Size interface for LinearMaps compatibility
Base.size(H::EffectiveHamiltonian) = (H.n, H.n)
Base.size(H::EffectiveHamiltonian, d::Int) = d <= 2 ? H.n : 1
Base.eltype(H::EffectiveHamiltonian{TL}) where {TL} = eltype(TL)

"""
    mul!(y, H::EffectiveHamiltonian, x)

Compute y = H * x = RHS⁻¹ * (LHS * x).
"""
function LinearAlgebra.mul!(y::AbstractVector, H::EffectiveHamiltonian, x::AbstractVector)
    # tmp = LHS * x
    mul!(H.tmp, H.LHS, x)
    # y = RHS⁻¹ * tmp
    ldiv!(y, H.RHS_factorization, H.tmp)
    return y
end

# Functor interface for Krylov.jl
(H::EffectiveHamiltonian)(x::AbstractVector) = H * x

# Matrix-vector multiplication
function Base.:*(H::EffectiveHamiltonian, x::AbstractVector)
    y = similar(x)
    mul!(y, H, x)
    return y
end

# ============================================================================
# NegatedOperator: A = -H for Krylov.jl interface
# ============================================================================

"""
    NegatedOperator{T}

Wrapper that negates an operator: A = -H.

This transforms the Green's function problem:
    G(z) = (z·I - H)⁻¹  →  solve (z·I + A)x = b where A = -H

Krylov.jl's `cg_lanczos_shift` solves `(A + σI)x = b`, so we use:
    A = -H, σ = z  →  (z·I - H)x = b

# Usage
```julia
H = EffectiveHamiltonian(LHS, RHS)
A = NegatedOperator(H)
# Now use A with Krylov.cg_lanczos_shift
```
"""
struct NegatedOperator{T}
    op::T
end

Base.size(A::NegatedOperator) = size(A.op)
Base.size(A::NegatedOperator, d::Int) = size(A.op, d)
Base.eltype(A::NegatedOperator) = eltype(A.op)

function LinearAlgebra.mul!(y::AbstractVector, A::NegatedOperator, x::AbstractVector)
    mul!(y, A.op, x)
    y .*= -1
    return y
end

(A::NegatedOperator)(x::AbstractVector) = A * x

function Base.:*(A::NegatedOperator, x::AbstractVector)
    y = A.op * x
    y .*= -1
    return y
end

# ============================================================================
# Matrix-Free Effective Hamiltonian: H = RHS⁻¹ * LHS
# ============================================================================

"""
    MatrixFreeEffectiveHamiltonian{D, W, T}

Matrix-free effective Hamiltonian operator H = RHS⁻¹ * LHS.

This is the fully matrix-free version of `EffectiveHamiltonian`, using O(N) memory
instead of O(N²) for dense matrices.

For photonic crystals:
- TE: RHS = μ, so RHS⁻¹ = 1/μ (element-wise in real space)
- TM: RHS = ε, so RHS⁻¹ = 1/ε (element-wise in real space)

For phononic crystals:
- SH/PSV: RHS = ρ, so RHS⁻¹ = 1/ρ (element-wise in real space)

# Grid Resolution Constraint

**Important**: The solver's grid resolution must be large enough to contain all
plane wave indices. The required resolution in each dimension is:

    resolution >= 2 * max_index + 1

where `max_index` is the maximum absolute value of plane wave indices (typically
equal to the `cutoff` parameter). For example, with `cutoff=7` (indices -7 to 7),
the resolution must be at least 15 in each dimension.

If the resolution is too small, Fourier coefficients outside the grid will be
silently dropped, leading to incorrect results (errors can be ~100%).

# Usage
```julia
op = MatrixFreeOperator(solver, k)
H = MatrixFreeEffectiveHamiltonian(op)
A = NegatedOperator(H)  # A = -H for RSCG
x, stats = rscg(A, b, shifts)
```

# RHS⁻¹ methods

- `:approximate`: Element-wise 1/ε in real space (fast, approximate for inhomogeneous media)
- `:cg`: Iterative CG to solve RHS * y = x (slower, exact)
"""
struct MatrixFreeEffectiveHamiltonian{
    D<:Dimension,W<:WaveType,T<:Complex,N,F,I,M<:RHSInvMethod
}
    op::MatrixFreeOperator{D,W,T,N,F,I}
    rhs_inv::Array{T}  # 1/material in real space (for ApproximateRHSInv)
    tmp_fourier::Vector{T}  # temporary for Fourier coefficients
    tmp_fourier2::Vector{T}  # temporary for CG
    tmp_grid::Array{T}  # temporary for real space
    n::Int  # problem dimension
    rhs_inv_method::M  # RHSInvMethod subtype (ApproximateRHSInv or CGRHSInv)
end

"""
    MatrixFreeEffectiveHamiltonian(op::MatrixFreeOperator, rhs_inv_method=ApproximateRHSInv())

Create a matrix-free effective Hamiltonian from a MatrixFreeOperator.

# Arguments
- `op::MatrixFreeOperator`: The matrix-free operator
- `rhs_inv_method::RHSInvMethod`: Method for RHS⁻¹ application
  - `ApproximateRHSInv()`: Element-wise 1/ε in real space (fast, but approximate for inhomogeneous media)
  - `CGRHSInv(; atol, rtol, maxiter)`: Iterative CG to solve RHS * y = x (slower, but exact)

# Examples
```julia
# Fast approximate method (default)
H = MatrixFreeEffectiveHamiltonian(op)
H = MatrixFreeEffectiveHamiltonian(op, ApproximateRHSInv())

# Accurate CG method for high-contrast media
H = MatrixFreeEffectiveHamiltonian(op, CGRHSInv())
H = MatrixFreeEffectiveHamiltonian(op, CGRHSInv(atol=1e-12, rtol=1e-12, maxiter=200))
```

# Backward Compatibility
Symbol-based `:approximate` and `:cg` are deprecated but still supported:
```julia
H = MatrixFreeEffectiveHamiltonian(op, :approximate)  # deprecated
H = MatrixFreeEffectiveHamiltonian(op, :cg)           # deprecated
```
"""
function MatrixFreeEffectiveHamiltonian(
    op::MatrixFreeOperator{D,W,T,N,F,I}, rhs_inv_method::RHSInvMethod=ApproximateRHSInv()
) where {D<:Dimension,W<:WaveType,T,N,F,I}
    solver = op.solver
    res = op.ctx.resolution
    nc = ncomponents(solver.wave)
    num_pw = solver.basis.num_pw
    n = num_pw * nc

    # Get RHS⁻¹ in real space (used for ApproximateRHSInv, and as initial guess for CGRHSInv)
    rhs_inv = _get_rhs_inv(solver.wave, solver.material_arrays, res, T)

    # Allocate temporaries
    tmp_fourier = zeros(T, n)
    tmp_fourier2 = zeros(T, n)
    tmp_grid = zeros(T, res)

    MatrixFreeEffectiveHamiltonian{D,W,T,N,F,I,typeof(rhs_inv_method)}(
        op, rhs_inv, tmp_fourier, tmp_fourier2, tmp_grid, n, rhs_inv_method
    )
end

# Backward compatibility: Symbol-based constructor (deprecated)
function MatrixFreeEffectiveHamiltonian(
    op::MatrixFreeOperator{D,W,T,N,F,I},
    rhs_inv_method::Symbol;
    cg_atol::Float64=1e-12,
    cg_rtol::Float64=1e-10,
    cg_maxiter::Int=100,
) where {D<:Dimension,W<:WaveType,T,N,F,I}
    Base.depwarn(
        "Symbol-based rhs_inv_method is deprecated. Use ApproximateRHSInv() or CGRHSInv() instead.",
        :MatrixFreeEffectiveHamiltonian,
    )
    method = if rhs_inv_method == :approximate
        ApproximateRHSInv()
    elseif rhs_inv_method == :cg
        CGRHSInv(; atol=cg_atol, rtol=cg_rtol, maxiter=cg_maxiter)
    else
        throw(ArgumentError("rhs_inv_method must be :approximate or :cg"))
    end
    MatrixFreeEffectiveHamiltonian(op, method)
end

# Helper to get 1/RHS material array
# Error fallback for unsupported wave types
function _get_rhs_inv(wave::WaveType, mats, res, ::Type{T}) where {T}
    throw(ArgumentError("Unsupported wave type for _get_rhs_inv: $(typeof(wave))"))
end

function _get_rhs_inv(::TEWave, mats, res, ::Type{T}) where {T}
    # TE: RHS = μ
    T.(1.0 ./ mats.μ)
end

function _get_rhs_inv(::TMWave, mats, res, ::Type{T}) where {T}
    # TM: RHS = ε
    T.(1.0 ./ mats.ε)
end

function _get_rhs_inv(::SHWave, mats, res, ::Type{T}) where {T}
    # SH: RHS = ρ
    T.(1.0 ./ mats.ρ)
end

function _get_rhs_inv(::PSVWave, mats, res, ::Type{T}) where {T}
    # PSV: RHS = ρ (same for both components)
    T.(1.0 ./ mats.ρ)
end

# Size interface
Base.size(H::MatrixFreeEffectiveHamiltonian) = (H.n, H.n)
Base.size(H::MatrixFreeEffectiveHamiltonian, d::Int) = d <= 2 ? H.n : 1
Base.eltype(::MatrixFreeEffectiveHamiltonian{D,W,T}) where {D,W,T} = T

"""
    mul!(y, H::MatrixFreeEffectiveHamiltonian, x)

Compute y = H * x = RHS⁻¹ * (LHS * x) in O(N log N) time.
"""
function LinearAlgebra.mul!(
    y::AbstractVector, H::MatrixFreeEffectiveHamiltonian{Dim2,W,T}, x::AbstractVector
) where {W<:WaveType,T}
    # Step 1: tmp = LHS * x
    apply_lhs!(H.tmp_fourier, H.op, x)

    # Step 2: y = RHS⁻¹ * tmp (element-wise division in real space)
    _apply_rhs_inv!(y, H, H.tmp_fourier)

    return y
end

# Helper for RHS⁻¹ application (uses trait-based + method-based dispatch)
function _apply_rhs_inv!(
    y::AbstractVector{T}, H::MatrixFreeEffectiveHamiltonian{Dim2,W,T}, x::AbstractVector{T}
) where {W<:WaveType,T}
    _apply_rhs_inv_impl!(wave_structure(W()), H.rhs_inv_method, y, H, x)
end

# Error fallback for unsupported combinations
function _apply_rhs_inv_impl!(trait, method, y, H, x)
    throw(
        ArgumentError(
            "Unsupported wave structure/method combination: $(typeof(trait)), $(typeof(method))",
        ),
    )
end

# ============================================================================
# ScalarWave2D implementations (TE, TM, SH)
# ============================================================================

# ApproximateRHSInv: element-wise 1/ε in real space (fast, approximate)
function _apply_rhs_inv_impl!(
    ::ScalarWave2D,
    ::ApproximateRHSInv,
    y::AbstractVector{T},
    H::MatrixFreeEffectiveHamiltonian{Dim2,W,T},
    x::AbstractVector{T},
) where {W<:WaveType,T}
    op = H.op
    basis = op.solver.basis
    res = resolution(op)

    # To real space
    fourier_to_grid!(H.tmp_grid, x, basis, res)

    # Multiply by 1/material
    H.tmp_grid .*= H.rhs_inv

    # Back to Fourier
    grid_to_fourier!(y, H.tmp_grid, basis, res)

    return y
end

# CGRHSInv: solve RHS * y = x iteratively (slower, exact)
function _apply_rhs_inv_impl!(
    ::ScalarWave2D,
    method::CGRHSInv,
    y::AbstractVector{T},
    H::MatrixFreeEffectiveHamiltonian{Dim2,W,T},
    x::AbstractVector{T},
) where {W<:WaveType,T}
    op = H.op

    # Solve RHS * y = x using CG
    # RHS is Hermitian positive definite, so CG is appropriate

    # Initial guess: use approximate solution as starting point
    _apply_rhs_inv_impl!(ScalarWave2D(), ApproximateRHSInv(), y, H, x)

    # r = x - RHS * y
    apply_rhs!(H.tmp_fourier2, op, y)
    r = x - H.tmp_fourier2

    # p = r
    p = copy(r)

    r_norm0 = norm(x)
    r_dot_r = real(dot(r, r))

    for iter in 1:method.maxiter
        # Check convergence
        r_norm = sqrt(r_dot_r)
        if r_norm < method.atol || r_norm / r_norm0 < method.rtol
            break
        end

        # Ap = RHS * p
        Ap = H.tmp_fourier2
        apply_rhs!(Ap, op, p)

        # α = (r, r) / (p, Ap)
        p_dot_Ap = real(dot(p, Ap))
        α = r_dot_r / p_dot_Ap

        # y = y + α * p
        y .+= α .* p

        # r = r - α * Ap
        r .-= α .* Ap

        r_dot_r_new = real(dot(r, r))

        # β = (r_new, r_new) / (r, r)
        β = r_dot_r_new / r_dot_r

        # p = r + β * p
        p .= r .+ β .* p

        r_dot_r = r_dot_r_new
    end

    return y
end

# ============================================================================
# VectorWave2D implementations (PSV - 2 components)
# ============================================================================

# ApproximateRHSInv for VectorWave2D
function _apply_rhs_inv_impl!(
    ::VectorWave2D,
    ::ApproximateRHSInv,
    y::AbstractVector{T},
    H::MatrixFreeEffectiveHamiltonian{Dim2,W,T},
    x::AbstractVector{T},
) where {W<:WaveType,T}
    op = H.op
    basis = op.solver.basis
    res = resolution(op)
    N = basis.num_pw

    # x and y components
    x_x = @view x[1:N]
    x_y = @view x[(N + 1):2N]
    y_x = @view y[1:N]
    y_y = @view y[(N + 1):2N]

    # x component
    fourier_to_grid!(H.tmp_grid, x_x, basis, res)
    H.tmp_grid .*= H.rhs_inv
    grid_to_fourier!(y_x, H.tmp_grid, basis, res)

    # y component
    fourier_to_grid!(H.tmp_grid, x_y, basis, res)
    H.tmp_grid .*= H.rhs_inv
    grid_to_fourier!(y_y, H.tmp_grid, basis, res)

    return y
end

# CGRHSInv for VectorWave2D
function _apply_rhs_inv_impl!(
    ::VectorWave2D,
    method::CGRHSInv,
    y::AbstractVector{T},
    H::MatrixFreeEffectiveHamiltonian{Dim2,W,T},
    x::AbstractVector{T},
) where {W<:WaveType,T}
    op = H.op

    # Solve RHS * y = x using CG (same structure as single component)
    # RHS is block-diagonal with ρ, so we can solve each component independently
    # But for generality, we solve the full system

    # Initial guess: use approximate solution as starting point
    _apply_rhs_inv_impl!(VectorWave2D(), ApproximateRHSInv(), y, H, x)

    # r = x - RHS * y
    apply_rhs!(H.tmp_fourier2, op, y)
    r = x - H.tmp_fourier2

    # p = r
    p = copy(r)

    r_norm0 = norm(x)
    r_dot_r = real(dot(r, r))

    for iter in 1:method.maxiter
        # Check convergence
        r_norm = sqrt(r_dot_r)
        if r_norm < method.atol || r_norm / r_norm0 < method.rtol
            break
        end

        # Ap = RHS * p
        Ap = H.tmp_fourier2
        apply_rhs!(Ap, op, p)

        # α = (r, r) / (p, Ap)
        p_dot_Ap = real(dot(p, Ap))
        α = r_dot_r / p_dot_Ap

        # y = y + α * p
        y .+= α .* p

        # r = r - α * Ap
        r .-= α .* Ap

        r_dot_r_new = real(dot(r, r))

        # β = (r_new, r_new) / (r, r)
        β = r_dot_r_new / r_dot_r

        # p = r + β * p
        p .= r .+ β .* p

        r_dot_r = r_dot_r_new
    end

    return y
end

# Functor interface for Krylov.jl
(H::MatrixFreeEffectiveHamiltonian)(x::AbstractVector) = H * x

function Base.:*(H::MatrixFreeEffectiveHamiltonian, x::AbstractVector)
    y = similar(x)
    mul!(y, H, x)
    return y
end

# ============================================================================
# 3D Matrix-Free Effective Hamiltonian
# ============================================================================

# Helper to get 1/RHS for 3D waves
function _get_rhs_inv(::FullVectorEM, mats, res, ::Type{T}) where {T}
    # FullVectorEM: RHS = μ (same for all 3 components)
    T.(1.0 ./ mats.μ)
end

function _get_rhs_inv(::FullElastic, mats, res, ::Type{T}) where {T}
    # FullElastic: RHS = ρ (same for all 3 components)
    T.(1.0 ./ mats.ρ)
end

# 3D mul! implementation
function LinearAlgebra.mul!(
    y::AbstractVector, H::MatrixFreeEffectiveHamiltonian{Dim3,W,T}, x::AbstractVector
) where {W<:WaveType,T}
    # Step 1: tmp = LHS * x
    apply_lhs!(H.tmp_fourier, H.op, x)

    # Step 2: y = RHS⁻¹ * tmp (element-wise division in real space)
    _apply_rhs_inv_3d!(y, H, H.tmp_fourier)

    return y
end

# Helper for RHS⁻¹ application for 3D waves (uses trait-based + method-based dispatch)
function _apply_rhs_inv_3d!(
    y::AbstractVector{T}, H::MatrixFreeEffectiveHamiltonian{Dim3,W,T}, x::AbstractVector{T}
) where {W<:WaveType,T}
    _apply_rhs_inv_3d_impl!(wave_structure(W()), H.rhs_inv_method, y, H, x)
end

# Error fallback for unsupported combinations (3D)
function _apply_rhs_inv_3d_impl!(trait, method, y, H, x)
    throw(
        ArgumentError(
            "Unsupported wave structure/method combination for 3D: $(typeof(trait)), $(typeof(method))",
        ),
    )
end

# ============================================================================
# VectorWave3D implementations (FullVectorEM, FullElastic - 3 components)
# ============================================================================

# ApproximateRHSInv for VectorWave3D
function _apply_rhs_inv_3d_impl!(
    ::VectorWave3D,
    ::ApproximateRHSInv,
    y::AbstractVector{T},
    H::MatrixFreeEffectiveHamiltonian{Dim3,W,T},
    x::AbstractVector{T},
) where {W<:WaveType,T}
    op = H.op
    basis = op.solver.basis
    res = resolution(op)
    N = basis.num_pw

    # 3 components: x, y, z
    x_x = @view x[1:N]
    x_y = @view x[(N + 1):2N]
    x_z = @view x[(2N + 1):3N]
    y_x = @view y[1:N]
    y_y = @view y[(N + 1):2N]
    y_z = @view y[(2N + 1):3N]

    # x component
    fourier_to_grid!(H.tmp_grid, x_x, basis, res)
    H.tmp_grid .*= H.rhs_inv
    grid_to_fourier!(y_x, H.tmp_grid, basis, res)

    # y component
    fourier_to_grid!(H.tmp_grid, x_y, basis, res)
    H.tmp_grid .*= H.rhs_inv
    grid_to_fourier!(y_y, H.tmp_grid, basis, res)

    # z component
    fourier_to_grid!(H.tmp_grid, x_z, basis, res)
    H.tmp_grid .*= H.rhs_inv
    grid_to_fourier!(y_z, H.tmp_grid, basis, res)

    return y
end

# CGRHSInv for VectorWave3D
function _apply_rhs_inv_3d_impl!(
    ::VectorWave3D,
    method::CGRHSInv,
    y::AbstractVector{T},
    H::MatrixFreeEffectiveHamiltonian{Dim3,W,T},
    x::AbstractVector{T},
) where {W<:WaveType,T}
    op = H.op

    # Solve RHS * y = x using CG
    # RHS is block-diagonal with μ (or ρ), so structure is similar to 2D

    # Initial guess: use approximate solution as starting point
    _apply_rhs_inv_3d_impl!(VectorWave3D(), ApproximateRHSInv(), y, H, x)

    # r = x - RHS * y
    apply_rhs!(H.tmp_fourier2, op, y)
    r = x - H.tmp_fourier2

    # p = r
    p = copy(r)

    r_norm0 = norm(x)
    r_dot_r = real(dot(r, r))

    for iter in 1:method.maxiter
        # Check convergence
        r_norm = sqrt(r_dot_r)
        if r_norm < method.atol || r_norm / r_norm0 < method.rtol
            break
        end

        # Ap = RHS * p
        Ap = H.tmp_fourier2
        apply_rhs!(Ap, op, p)

        # α = (r, r) / (p, Ap)
        p_dot_Ap = real(dot(p, Ap))
        α = r_dot_r / p_dot_Ap

        # y = y + α * p
        y .+= α .* p

        # r = r - α * Ap
        r .-= α .* Ap

        r_dot_r_new = real(dot(r, r))

        # β = (r_new, r_new) / (r, r)
        β = r_dot_r_new / r_dot_r

        # p = r + β * p
        p .= r .+ β .* p

        r_dot_r = r_dot_r_new
    end

    return y
end

# ============================================================================
# Unified Green's Function API
# ============================================================================

"""
    GFMethod

Abstract type for Green's function computation methods.
Used by `compute_greens_function`, `compute_dos`, and `compute_ldos`.
"""
abstract type GFMethod end

"""
    DirectGF()

Dense direct solve method (LU factorization).
Most accurate but O(N²) memory and O(N³) per frequency.
"""
struct DirectGF <: GFMethod end

"""
    RSKGF(; atol=1e-10, rtol=1e-10, itmax=0, verbose=0)

Dense RSK method using ReducedShiftedKrylov.jl.
O(N²) memory, efficient for many frequencies.
"""
struct RSKGF <: GFMethod
    atol::Float64
    rtol::Float64
    itmax::Int
    verbose::Int
end
function RSKGF(; atol::Float64=1e-10, rtol::Float64=1e-10, itmax::Int=0, verbose::Int=0)
    RSKGF(atol, rtol, itmax, verbose)
end

"""
    MatrixFreeGF(; rhs_inv_method=:approximate, atol=1e-10, rtol=1e-10, itmax=0, verbose=0)

Matrix-free RSCG method. O(N) memory, O(N log N) per iteration.

# Arguments
- `rhs_inv_method`: Method for RHS⁻¹ application
  - `:approximate`: Element-wise 1/ε (fast, approximate for inhomogeneous media)
  - `:cg`: Iterative CG (slower, exact)

# Recommended Usage
- Peak/resonance detection: accurate even without full RSCG convergence
- Large-scale calculations where Dense methods are memory-prohibitive
- For accurate absolute values, use `DirectGF()` instead
"""
struct MatrixFreeGF{M<:RHSInvMethod} <: GFMethod
    rhs_inv_method::M
    atol::Float64
    rtol::Float64
    itmax::Int
    verbose::Int
end

# Unified constructor - accepts both RHSInvMethod and Symbol (deprecated)
function MatrixFreeGF(;
    rhs_inv_method::Union{RHSInvMethod,Symbol}=ApproximateRHSInv(),
    atol::Float64=1e-10,
    rtol::Float64=1e-10,
    itmax::Int=0,
    verbose::Int=0,
)
    method = _convert_rhs_inv_method(rhs_inv_method)
    MatrixFreeGF{typeof(method)}(method, atol, rtol, itmax, verbose)
end

# ============================================================================
# Unified compute_greens_function
# ============================================================================

"""
    compute_greens_function(solver, k, ω_values, source, method::GFMethod; η=1e-3)

Compute Green's function G(ω) * source using the specified method.

# Methods
- `DirectGF()`: Dense direct solve (most accurate, O(N²) memory)
- `RSKGF()`: Dense RSK using ReducedShiftedKrylov.jl
- `MatrixFreeGF()`: Matrix-free RSCG (O(N) memory, best for large systems)

# Examples
```julia
# Direct solve (default)
G = compute_greens_function(solver, k, ω_values, source; η=1e-2)

# Matrix-free for large systems
G = compute_greens_function(solver, k, ω_values, source, MatrixFreeGF(); η=1e-2)
```
"""
function compute_greens_function(
    solver::Solver,
    k,
    ω_values::AbstractVector{<:Real},
    source::AbstractVector,
    method::DirectGF;
    η::Real=1e-3,
)
    N = solver.basis.num_pw
    T = ComplexF64

    LHS, RHS = build_matrices(solver, k)
    RHS_source = RHS * source

    results = Vector{Vector{T}}(undef, length(ω_values))

    for (i, ω) in enumerate(ω_values)
        z = ω^2 + im*η
        A = z * RHS - LHS
        results[i] = A \ RHS_source
    end

    return results
end

# Generic fallback for GFMethod - extension will provide specific implementations
# Note: Do NOT add specific fallbacks for RSKGF/MatrixFreeGF here as they would
# conflict with package extension methods (Julia does not allow method overwriting)
function compute_greens_function(
    solver::Solver,
    k,
    ω_values::AbstractVector{<:Real},
    source::AbstractVector,
    method::GFMethod;
    η::Real=1e-3,
)
    method_name = typeof(method).name.name
    if method_name in (:RSKGF, :MatrixFreeGF)
        throw(
            ArgumentError(
                "$method_name requires ReducedShiftedKrylov.jl. " *
                "Add `using ReducedShiftedKrylov` to enable this feature.",
            ),
        )
    else
        throw(
            ArgumentError(
                "Unsupported GFMethod: $(typeof(method)). " *
                "Available methods: DirectGF(), RSKGF(), MatrixFreeGF()",
            ),
        )
    end
end

# ============================================================================
# Unified compute_dos
# ============================================================================

"""
    compute_dos(solver, ω_values, k_points, method::GFMethod; η=1e-3, n_random=10)

Compute density of states (DOS) using the specified method.

# Methods
- `DirectGF()`: Dense direct solve (most accurate)
- `RSKGF()`: Dense RSK with stochastic trace estimation
- `MatrixFreeGF()`: Matrix-free RSCG with stochastic trace estimation

# Examples
```julia
# Direct solve (default)
dos = compute_dos(solver, ω_values, k_points; η=1e-2)

# Matrix-free for large systems
dos = compute_dos(solver, ω_values, k_points, MatrixFreeGF(); η=1e-2)
```
"""
function compute_dos(
    solver::Solver{Dim2},
    ω_values::AbstractVector{<:Real},
    k_points::AbstractVector,
    method::DirectGF;
    η::Real=1e-3,
)
    N = solver.basis.num_pw
    dos = zeros(length(ω_values))

    for k in k_points
        k_vec = _to_kvec(k)
        LHS, RHS = build_matrices(solver, k_vec)

        for (iω, ω) in enumerate(ω_values)
            z = ω^2 + im*η
            A = z * RHS - LHS
            A_factored = lu(A)
            trace_G = zero(ComplexF64)

            for n in 1:N
                e_n = zeros(ComplexF64, N)
                e_n[n] = 1.0
                b = RHS * e_n
                x = A_factored \ b
                trace_G += x[n]
            end

            dos[iω] += -imag(trace_G) / π
        end
    end

    dos ./= length(k_points)
    return dos
end

# Generic fallback for GFMethod
function compute_dos(
    solver::Solver,
    ω_values::AbstractVector{<:Real},
    k_points::AbstractVector,
    method::GFMethod;
    η::Real=1e-3,
    n_random::Int=10,
)
    method_name = typeof(method).name.name
    if method_name in (:RSKGF, :MatrixFreeGF)
        throw(
            ArgumentError(
                "$method_name requires ReducedShiftedKrylov.jl. " *
                "Add `using ReducedShiftedKrylov` to enable this feature.",
            ),
        )
    else
        throw(
            ArgumentError(
                "Unsupported GFMethod: $(typeof(method)). " *
                "Available methods: DirectGF(), RSKGF(), MatrixFreeGF()",
            ),
        )
    end
end

# ============================================================================
# Unified compute_ldos
# ============================================================================

"""
    compute_ldos(solver, position, ω_values, k_points, method::GFMethod; η=1e-3)

Compute local density of states (LDOS) using the specified method.

# Methods
- `DirectGF()`: Dense direct solve (most accurate, O(N²) memory)
- `RSKGF()`: Dense RSK using ReducedShiftedKrylov.jl
- `MatrixFreeGF()`: Matrix-free RSCG (O(N) memory, best for large systems)

# Examples
```julia
# Direct solve (default when method not specified)
ldos = compute_ldos(solver, pos, ω_values, k_points; η=1e-2)

# Matrix-free for large systems
ldos = compute_ldos(solver, pos, ω_values, k_points, MatrixFreeGF(); η=1e-2)

# Matrix-free with exact RHS⁻¹
ldos = compute_ldos(solver, pos, ω_values, k_points, MatrixFreeGF(rhs_inv_method=:cg); η=1e-2)

# Dense RSK
ldos = compute_ldos(solver, pos, ω_values, k_points, RSKGF(); η=1e-2)
```

# See Also
- `compute_dos`: Compute density of states (trace of Green's function)
- `compute_greens_function`: Compute Green's function directly
"""
function compute_ldos(
    solver::Solver{Dim2},
    position::AbstractVector{<:Real},
    ω_values::AbstractVector{<:Real},
    k_points::AbstractVector,
    method::DirectGF;
    η::Real=1e-3,
)
    N = solver.basis.num_pw
    ldos = zeros(length(ω_values))

    r = position
    source = [exp(-im * dot(G, r)) for G in solver.basis.G]
    source ./= N

    for k in k_points
        k_vec = _to_kvec(k)
        G_values = compute_greens_function(solver, k_vec, ω_values, source; η=η)

        for (iω, G) in enumerate(G_values)
            G_rr = dot(conj.(source), G)
            ldos[iω] += -imag(G_rr) / π
        end
    end

    ldos ./= length(k_points)
    return ldos
end

# Generic fallback for GFMethod
function compute_ldos(
    solver::Solver,
    position::AbstractVector{<:Real},
    ω_values::AbstractVector{<:Real},
    k_points::AbstractVector,
    method::GFMethod;
    η::Real=1e-3,
)
    method_name = typeof(method).name.name
    if method_name in (:RSKGF, :MatrixFreeGF)
        throw(
            ArgumentError(
                "$method_name requires ReducedShiftedKrylov.jl. " *
                "Add `using ReducedShiftedKrylov` to enable this feature.",
            ),
        )
    else
        throw(
            ArgumentError(
                "Unsupported GFMethod: $(typeof(method)). " *
                "Available methods: DirectGF(), RSKGF(), MatrixFreeGF()",
            ),
        )
    end
end
