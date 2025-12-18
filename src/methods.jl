# Last-Modified: 2025-12-12T22:18:28+09:00

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
end

"""
    KrylovKitMethod(; tol=1e-8, maxiter=300, krylovdim=30, verbosity=0, shift=0.0)

Create a KrylovKit-based iterative eigenvalue solver.

# Keyword Arguments
- `tol`: Convergence tolerance (default: 1e-8)
- `maxiter`: Maximum number of iterations (default: 300)
- `krylovdim`: Dimension of Krylov subspace (default: 30)
- `verbosity`: Output level, 0=silent, 1=warnings, 2=info (default: 0)
- `shift`: Spectral shift for shift-and-invert transformation (default: 0.0)
  - Use `shift=0.01` for 3D photonic crystals to skip spurious longitudinal modes
  - Use larger values to target specific frequency ranges

# Example
```julia
# Standard iterative solver
solver = Solver(TEWave(), geo, (64, 64), KrylovKitMethod(); cutoff=10)

# 3D with shift-and-invert to skip longitudinal modes
solver = Solver(FullVectorEM(), geo, (16, 16, 16), KrylovKitMethod(shift=0.01); cutoff=3)
```
"""
function KrylovKitMethod(; tol::Real=1e-8, maxiter::Int=300, krylovdim::Int=30,
                          verbosity::Int=0, shift::Real=0.0)
    KrylovKitMethod(Float64(tol), maxiter, krylovdim, verbosity, Float64(shift))
end

"""
    LOBPCGMethod <: IterativeMethod

Locally Optimal Block Preconditioned Conjugate Gradient method.
Based on Knyazev (2001), SIAM J. Sci. Comput. Vol.23, No.2, pp.517-541.

Solves the symmetric generalized eigenvalue problem `A x = λ B x` where
both A and B are Hermitian and B is positive definite. This method is
particularly effective for:
- Large-scale problems with many eigenvalues
- Problems where good preconditioners are available
- Phononic crystal calculations where eigenvalue scaling may be an issue

# Fields
- `tol::Float64`: Convergence tolerance (default: 1e-8)
- `maxiter::Int`: Maximum iterations (default: 200)
- `shift::Float64`: Spectral shift for shift-and-invert transformation (default: 0.0)

# Notes
- Requires symmetric/Hermitian matrices A and B
- B must be positive definite
- Works directly with dense matrices (no matrix-free support yet)
- Block method: computes multiple eigenvalues simultaneously
- When shift > 0, uses shift-and-invert to skip eigenvalues near zero
  (useful for 3D H-field formulation where spurious modes exist at λ ≈ 0)

# Example
```julia
# For phononic calculations
solver = Solver(SHWave(), geo, (64, 64), LOBPCGMethod())

# With tighter tolerance
method = LOBPCGMethod(tol=1e-10, maxiter=500)

# For 3D FullVectorEM (skip spurious modes at λ ≈ 0)
method = LOBPCGMethod(shift=0.01)
```

See also: [`KrylovKitMethod`](@ref), [`DenseMethod`](@ref)
"""
struct LOBPCGMethod <: IterativeMethod
    tol::Float64
    maxiter::Int
    shift::Float64
end

"""
    LOBPCGMethod(; tol=1e-8, maxiter=200, shift=0.0)

Create a LOBPCG-based iterative eigenvalue solver.

# Keyword Arguments
- `tol`: Convergence tolerance (default: 1e-8)
- `maxiter`: Maximum number of iterations (default: 200)
- `shift`: Spectral shift σ for shift-and-invert (default: 0.0)
  - When shift > 0, transforms `A x = λ B x` to `(A - σB)⁻¹ B x = μ x`
  - Only eigenvalues λ > σ are returned
  - Useful for 3D FullVectorEM to skip spurious longitudinal modes at λ ≈ 0

# Example
```julia
solver = Solver(PSVWave(), geo, (64, 64), LOBPCGMethod())

# For 3D with shift-and-invert
solver = Solver(FullVectorEM(), geo, (16, 16, 16), LOBPCGMethod(shift=0.01); cutoff=3)
```
"""
function LOBPCGMethod(; tol::Real=1e-8, maxiter::Int=200, shift::Real=0.0)
    LOBPCGMethod(Float64(tol), maxiter, Float64(shift))
end

# Future extensions:
# struct SeededRSCG <: RSCGMethod end   # Uses previous solution as seed
# struct BlockRSCG <: RSCGMethod end    # Block version for multiple sources
