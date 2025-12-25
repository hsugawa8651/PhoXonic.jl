# Last-Modified: 2025-12-25T12:30:00+09:00

#=
Solver methods for eigenvalue problems.
Defines the type hierarchy for different solution approaches.
=#

"""
    SolverMethod

Abstract type for solver methods.
"""
abstract type SolverMethod end

"""
    IterativeMethod <: SolverMethod

Abstract type for iterative solver methods.
"""
abstract type IterativeMethod <: SolverMethod end

"""
    RSCGMethod <: IterativeMethod

Abstract type for Reduced Shifted Conjugate Gradient methods.
Used for Green's function calculations.
"""
abstract type RSCGMethod <: IterativeMethod end

# ============================================================================
# Concrete types
# ============================================================================

"""
    DenseMethod <: SolverMethod

Dense matrix eigenvalue solver using LAPACK's `eigen` function.

This is the default method, suitable for small to medium-sized systems.
It builds the full N×N matrices and computes all eigenvalues/eigenvectors.

# Fields
- `shift::Float64`: Minimum eigenvalue cutoff (default: 0.0)
  - Eigenvalues with ω² < shift are filtered out after solving
  - Use `shift > 0` for 3D FullVectorEM to skip spurious longitudinal modes (λ ≈ 0)
  - Recommended: `shift = 0.01` for 3D photonic crystals
  - Note: Unlike iterative methods, this is post-hoc filtering, not shift-and-invert

# Complexity
- Memory: O(N²)
- Time: O(N³)

where N = num_plane_waves × ncomponents(wave)

# Recommended For
- 1D systems (any resolution)
- 2D systems with N < ~1000
- Development and debugging (exact results)

# Example
```julia
# Default (implicit)
solver = Solver(TEWave(), geo, (64, 64))

# Explicit
solver = Solver(TEWave(), geo, (64, 64), DenseMethod())

# 3D with shift to skip spurious modes
solver = Solver(FullVectorEM(), geo, (12, 12, 12), DenseMethod(shift=0.01))
```

See also: [`KrylovKitMethod`](@ref), [`LOBPCGMethod`](@ref)
"""
struct DenseMethod <: SolverMethod
    shift::Float64
end

"""
    DenseMethod(; shift=0.0)

Create a dense matrix eigenvalue solver.

# Keyword Arguments
- `shift`: Minimum eigenvalue cutoff (default: 0.0)
  - Eigenvalues with ω² < shift are filtered out
  - Use for 3D FullVectorEM to skip spurious longitudinal modes
  - API is consistent with KrylovKitMethod and LOBPCGMethod

# Example
```julia
# Standard usage (2D)
solver = Solver(TEWave(), geo, (64, 64), DenseMethod())

# 3D with shift to skip spurious modes
solver = Solver(FullVectorEM(), geo, (12, 12, 12), DenseMethod(shift=0.01))
```
"""
DenseMethod(; shift::Real=0.0) = DenseMethod(Float64(shift))

"""
    BasicRSCG <: RSCGMethod

Basic RSCG method for Green's function computation.
Suitable for DOS/LDOS calculations in large systems.
"""
struct BasicRSCG <: RSCGMethod end

"""
    KrylovKitMethod <: IterativeMethod

Iterative eigenvalue solver using KrylovKit.jl.
Suitable for large systems where dense matrix storage is impractical.

Uses shift-and-invert spectral transformation when `shift > 0`:
- Original problem: `A x = λ B x`
- Transformed: `(A - σB)⁻¹ B x = μ x` where `μ = 1/(λ - σ)`

This is particularly useful for 3D calculations where the H-field formulation
produces N spurious longitudinal modes with λ ≈ 0. Setting `shift > 0` (e.g., 0.01)
effectively filters out these unphysical modes and focuses on the transverse
(physical) modes.

# Fields
- `tol::Float64`: Convergence tolerance (default: 1e-8)
- `maxiter::Int`: Maximum iterations (default: 300)
- `krylovdim::Int`: Krylov subspace dimension (default: 30)
- `verbosity::Int`: Output verbosity level (0=silent, 1=warn, 2=info)
- `shift::Float64`: Spectral shift σ for shift-and-invert (default: 0.0)
  - `shift = 0`: Standard generalized eigenvalue problem (no transformation)
  - `shift > 0`: Shift-and-invert, finds eigenvalues closest to σ
  - Recommended: `shift = 0.01` for 3D photonic crystals to skip longitudinal modes
- `matrix_free::Bool`: Use matrix-free operators for shift-and-invert (default: false)
  - `false`: Build dense matrices (faster for N < 2000)
  - `true`: Use O(N) memory with iterative inner solver (for large 3D problems)

# Example
```julia
# 2D calculation (no shift needed)
method = KrylovKitMethod()

# 3D calculation (use shift to skip longitudinal modes)
method = KrylovKitMethod(shift=0.01)

# Target specific frequency range
method = KrylovKitMethod(shift=1.5)  # Find modes near ω² = 1.5
```

See also: [`DenseMethod`](@ref), [`Solver`](@ref)
"""
struct KrylovKitMethod <: IterativeMethod
    tol::Float64
    maxiter::Int
    krylovdim::Int
    verbosity::Int
    shift::Float64
    matrix_free::Bool
end

"""
    KrylovKitMethod(; tol=1e-8, maxiter=300, krylovdim=30, verbosity=0, shift=0.0, matrix_free=false)

Create a KrylovKit-based iterative eigenvalue solver.

# Keyword Arguments
- `tol`: Convergence tolerance (default: 1e-8)
- `maxiter`: Maximum number of iterations (default: 300)
- `krylovdim`: Dimension of Krylov subspace (default: 30)
- `verbosity`: Output level, 0=silent, 1=warnings, 2=info (default: 0)
- `shift`: Spectral shift for shift-and-invert transformation (default: 0.0)
  - Use `shift=0.01` for 3D photonic crystals to skip spurious longitudinal modes
  - Use larger values to target specific frequency ranges
- `matrix_free`: Use matrix-free operators for shift-and-invert (default: false)
  - `false`: Build dense matrices (faster for N < 2000, O(N²) memory)
  - `true`: Use iterative inner solver (O(N) memory, for large 3D problems)

# Example
```julia
# Standard iterative solver
solver = Solver(TEWave(), geo, (64, 64), KrylovKitMethod(); cutoff=10)

# 3D with shift-and-invert to skip longitudinal modes
solver = Solver(FullVectorEM(), geo, (16, 16, 16), KrylovKitMethod(shift=0.01); cutoff=3)

# Large 3D with matrix-free shift-and-invert (memory-efficient)
solver = Solver(FullVectorEM(), geo, (32, 32, 32), KrylovKitMethod(shift=0.01, matrix_free=true); cutoff=5)
```
"""
function KrylovKitMethod(;
    tol::Real=1e-8,
    maxiter::Int=300,
    krylovdim::Int=30,
    verbosity::Int=0,
    shift::Real=0.0,
    matrix_free::Bool=false,
)
    KrylovKitMethod(
        Float64(tol), maxiter, krylovdim, verbosity, Float64(shift), matrix_free
    )
end

"""
    LOBPCGMethod <: IterativeMethod

Locally Optimal Block Preconditioned Conjugate Gradient method.
Based on Knyazev (2001), SIAM J. Sci. Comput. Vol.23, No.2, pp.517-541.

Solves the symmetric generalized eigenvalue problem `A x = λ B x` where
both A and B are Hermitian and B is positive definite. This method is
particularly effective for:
- Large-scale band structure calculations with `compute_bands`
- Phononic crystal calculations (steel/epoxy, silicon/air, etc.)

With warm start enabled (default), LOBPCG reuses eigenvectors from previous
k-points as initial guesses, achieving significant speedup while maintaining
accuracy. The first k-point uses Dense (if first_dense=true) to provide
accurate initial eigenvectors.

# Fields
- `tol::Float64`: Convergence tolerance (default: 1e-3)
- `maxiter::Int`: Maximum iterations (default: 100)
- `shift::Float64`: Spectral shift for shift-and-invert transformation (default: 0.0)
- `warm_start::Bool`: Use previous eigenvectors as initial guess (default: true)
- `scale::Bool`: Scale matrix A by max|A| (default: false, scaling can hurt phononic problems)
- `first_dense::Bool`: Solve first k-point with Dense for accurate warm start (default: true)
- `preconditioner`: Preconditioner type (default: :none)

# Notes
- Requires symmetric/Hermitian matrices A and B
- B must be positive definite
- Works directly with dense matrices (no matrix-free support yet)
- Block method: computes multiple eigenvalues simultaneously
- For phononic problems, leave scale=false and preconditioner=:none (defaults)
- When shift > 0, uses shift-and-invert to skip eigenvalues near zero
  (useful for 3D H-field formulation where spurious modes exist at λ ≈ 0)
- **Note**: With shift > 0, the method falls back to dense `eigen` solver internally
  because the shifted matrix (A - σB) is not positive definite as required by LOBPCG.

# Example
```julia
# Default with warm start (recommended for band structure)
solver = Solver(PSVWave(), geo, (64, 64), LOBPCGMethod(); cutoff=20)

# For 3D FullVectorEM (skip spurious modes at λ ≈ 0)
method = LOBPCGMethod(shift=0.01)
```

See also: [`KrylovKitMethod`](@ref), [`DenseMethod`](@ref), [`matrix_dimension`](@ref)
"""
struct LOBPCGMethod <: IterativeMethod
    tol::Float64
    maxiter::Int
    shift::Float64
    warm_start::Bool
    scale::Bool
    first_dense::Bool
    preconditioner::Any  # Symbol or custom preconditioner object
end

"""
    LOBPCGMethod(; tol=1e-3, maxiter=100, shift=0.0, warm_start=true, scale=false, first_dense=true, preconditioner=:none)

Create a LOBPCG-based iterative eigenvalue solver.

# Keyword Arguments
- `tol`: Convergence tolerance (default: 1e-3)
- `maxiter`: Maximum number of iterations (default: 100)
- `shift`: Spectral shift σ for shift-and-invert (default: 0.0)
  - When shift > 0, transforms `A x = λ B x` to `(A - σB)⁻¹ B x = μ x`
  - Only eigenvalues λ > σ are returned
  - Useful for 3D FullVectorEM to skip spurious longitudinal modes at λ ≈ 0
- `warm_start`: Use previous k-point's eigenvectors as initial guess (default: true)
- `scale`: Scale matrix A by max|A| (default: false)
  - For phononic problems, scaling can hurt convergence; leave as false
  - For photonic problems with normalized units, scaling may help
- `first_dense`: Solve first k-point with Dense for accurate initial eigenvectors (default: true)
- `preconditioner`: Preconditioner type (default: :none)
  - `:none`: No preconditioner (recommended for phononic problems)
  - `:diagonal`: Diagonal preconditioner using diag(A)^{-1}
  - Custom object: Any object implementing `ldiv!(y, P, x)`
  - Note: diagonal preconditioner may cause convergence issues for phononic problems

# Example
```julia
# Default (warm start enabled, recommended for band structure)
solver = Solver(PSVWave(), geo, (64, 64), LOBPCGMethod(); cutoff=20)

# For 3D with shift-and-invert
solver = Solver(FullVectorEM(), geo, (16, 16, 16), LOBPCGMethod(shift=0.01); cutoff=3)
```
"""
function LOBPCGMethod(;
    tol::Real=1e-3,
    maxiter::Int=100,
    shift::Real=0.0,
    warm_start::Bool=true,
    scale::Bool=false,
    first_dense::Bool=true,
    preconditioner=:none,
)
    LOBPCGMethod(
        Float64(tol),
        maxiter,
        Float64(shift),
        warm_start,
        scale,
        first_dense,
        preconditioner,
    )
end

# Future extensions:
# struct SeededRSCG <: RSCGMethod end   # Uses previous solution as seed
# struct BlockRSCG <: RSCGMethod end    # Block version for multiple sources
