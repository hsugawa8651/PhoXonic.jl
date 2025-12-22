# Last-Modified: 2025-12-22T23:00:00+09:00

#=
Eigenvalue problem solver for photonic and phononic crystals.
=#

"""
    AbstractSolver

Abstract base type for all solvers in PhoXonic.jl.

Subtypes:
- `Solver`: Plane wave expansion (PWE) solver

This abstract type allows for future extensions (e.g., GME solver).
"""
abstract type AbstractSolver end

"""
    Solver{D<:Dimension, W<:WaveType, M<:SolverMethod} <: AbstractSolver

Solver for computing band structures of photonic/phononic crystals using
the plane wave expansion (PWE) method.

# Fields
- `wave`: Wave type (TE, TM, SH, etc.)
- `geometry`: Crystal geometry
- `basis`: Plane wave basis
- `resolution`: Grid resolution for discretization
- `material_arrays`: Pre-computed material arrays on the grid
- `method`: Solver method (DenseMethod(), BasicRSCG(), etc.)
"""
struct Solver{D<:Dimension,W<:WaveType,M<:SolverMethod} <: AbstractSolver
    wave::W
    geometry::Geometry{D}
    basis::PlaneWaveBasis{D}
    resolution::NTuple{N,Int} where {N}
    material_arrays::NamedTuple
    method::M
end

"""
    Solver(wave, geometry, resolution, [method]; cutoff=7, discretization=SimpleGrid())

Create a solver for the given wave type and geometry.

# Arguments
- `wave`: Wave type (TEWave(), TMWave(), SHWave(), etc.)
- `geometry`: Crystal geometry
- `resolution`: Grid resolution (e.g., (64, 64) for 2D)
- `method`: Solver method (default: DenseMethod())

# Keyword Arguments
- `cutoff`: Plane wave cutoff (default: 7)
- `discretization`: Discretization method (SimpleGrid() or SubpixelAverage())

# Examples
```julia
# Default (Dense)
solver = Solver(TEWave(), geo, (64, 64))

# With explicit method
solver = Solver(TEWave(), geo, (64, 64), KrylovKitMethod())

# With cutoff
solver = Solver(TEWave(), geo, (64, 64), LOBPCGMethod(); cutoff=10)
```
"""
function Solver(
    wave::W,
    geometry::Geometry{D},
    resolution::NTuple{N,Int},
    method::M=DenseMethod();
    cutoff::Int=7,
    discretization::DiscretizationMethod=SimpleGrid(),
) where {D<:Dimension,W<:WaveType,N,M<:SolverMethod}
    basis = PlaneWaveBasis(geometry.lattice, cutoff)
    material_arrays = prepare_materials(wave, geometry, resolution, discretization)
    Solver{D,W,M}(wave, geometry, basis, resolution, material_arrays, method)
end

# Convenience constructor for 1D: accept single Int instead of tuple
function Solver(
    wave::W,
    geometry::Geometry{Dim1},
    resolution::Int,
    method::M=DenseMethod();
    cutoff::Int=7,
    discretization::DiscretizationMethod=SimpleGrid(),
) where {W<:WaveType,M<:SolverMethod}
    Solver(
        wave, geometry, (resolution,), method; cutoff=cutoff, discretization=discretization
    )
end

# ============================================================================
# Matrix dimension helper
# ============================================================================

"""
    matrix_dimension(solver::Solver) -> Int

Return the dimension of the eigenvalue problem matrix.

This is the size of the matrices A and B in the generalized eigenvalue problem
`A x = λ B x`. Use this to determine the size of initial vectors for iterative
solvers like LOBPCG.

# Returns
- `num_pw` for scalar waves (TE, TM, SH, 1D waves)
- `2 * num_pw` for 2D vector waves (P-SV)
- `3 * num_pw` for 3D vector waves (FullVectorEM, FullElastic)

# Example
```julia
solver = Solver(PSVWave(), geo, (64, 64); cutoff=20)
dim = matrix_dimension(solver)  # 2 * num_pw

# Create initial vectors for LOBPCG
X0 = randn(ComplexF64, dim, 20)
X0, _ = qr(X0)
X0 = Matrix(X0)
```

"""
function matrix_dimension(solver::Solver)
    ncomponents(solver.wave) * solver.basis.num_pw
end

# ============================================================================
# Bands selection helpers (multiple dispatch instead of isa checks)
# ============================================================================

"""
    _nev(bands, dim; default_max=20) -> Int

Compute number of eigenvalues to request from iterative solver.
"""
_nev(::Colon, dim; default_max::Int=20) = min(dim, default_max)
_nev(bands, dim; default_max::Int=20) = maximum(bands)

"""
    _select_bands(frequencies, eigenvectors, bands)

Select requested bands from computed results.
"""
_select_bands(frequencies, eigenvectors, ::Colon) = (frequencies, eigenvectors)
function _select_bands(frequencies, eigenvectors, bands)
    (frequencies[bands], eigenvectors[:, bands])
end

# Variant that returns first n bands (for KrylovKit where we get all eigenvalues sorted)
_select_first_bands(frequencies, eigenvectors, ::Colon) = (frequencies, eigenvectors)
function _select_first_bands(frequencies, eigenvectors, bands)
    nbands = min(length(frequencies), maximum(bands))
    (frequencies[1:nbands], eigenvectors[:, 1:nbands])
end

# Variant that filters available bands with warning support (for LOBPCG)
function _select_available_bands(frequencies, eigenvectors, ::Colon, _)
    (frequencies, eigenvectors, false)
end
function _select_available_bands(frequencies, eigenvectors, bands, nev_available)
    actual_bands = filter(b -> b <= nev_available, bands)
    warn_missing = length(actual_bands) < length(bands)
    (frequencies[actual_bands], eigenvectors[:, actual_bands], warn_missing)
end

# ============================================================================
# Material preparation
# ============================================================================

"""
    prepare_materials(wave::WaveType, geometry::Geometry, resolution, method)

Prepare discretized material arrays for the solver.

# Arguments
- `wave`: Wave type
- `geometry`: Crystal geometry
- `resolution`: Grid resolution
- `method`: Discretization method (SimpleGrid or SubpixelAverage)
"""
# Error fallback for unsupported wave type / dimension combinations
function prepare_materials(
    wave::WaveType,
    geometry::Geometry,
    resolution,
    method::DiscretizationMethod=SimpleGrid(),
)
    throw(
        ArgumentError(
            "Unsupported wave type / dimension combination: $(typeof(wave)) with $(typeof(geometry).parameters[1])",
        ),
    )
end

function prepare_materials(
    ::TEWave,
    geometry::Geometry{Dim2},
    resolution,
    method::DiscretizationMethod=SimpleGrid(),
)
    ε = discretize(geometry, resolution, :ε, method)
    μ = discretize(geometry, resolution, :μ, method)
    # For TE: use properly averaged ε⁻¹ (important for SubpixelAverage)
    ε_inv = discretize(geometry, resolution, :ε_inv, method)
    (ε=ε, μ=μ, ε_inv=ε_inv)
end

function prepare_materials(
    ::TMWave,
    geometry::Geometry{Dim2},
    resolution,
    method::DiscretizationMethod=SimpleGrid(),
)
    ε = discretize(geometry, resolution, :ε, method)
    μ = discretize(geometry, resolution, :μ, method)
    # For TM: use properly averaged μ⁻¹
    μ_inv = discretize(geometry, resolution, :μ_inv, method)
    (ε=ε, μ=μ, μ_inv=μ_inv)
end

function prepare_materials(
    ::SHWave,
    geometry::Geometry{Dim2},
    resolution,
    method::DiscretizationMethod=SimpleGrid(),
)
    ρ = discretize(geometry, resolution, :ρ, method)
    C44 = discretize(geometry, resolution, :C44, method)
    (ρ=ρ, C44=C44)
end

function prepare_materials(
    ::PSVWave,
    geometry::Geometry{Dim2},
    resolution,
    method::DiscretizationMethod=SimpleGrid(),
)
    ρ = discretize(geometry, resolution, :ρ, method)
    C11 = discretize(geometry, resolution, :C11, method)
    C12 = discretize(geometry, resolution, :C12, method)
    C44 = discretize(geometry, resolution, :C44, method)
    (ρ=ρ, C11=C11, C12=C12, C44=C44)
end

# 1D Photonic
function prepare_materials(
    ::Photonic1D,
    geometry::Geometry{Dim1},
    resolution,
    method::DiscretizationMethod=SimpleGrid(),
)
    ε = discretize(geometry, resolution, :ε, method)
    μ = discretize(geometry, resolution, :μ, method)
    ε_inv = discretize(geometry, resolution, :ε_inv, method)
    (ε=ε, μ=μ, ε_inv=ε_inv)
end

# 1D Longitudinal elastic
function prepare_materials(
    ::Longitudinal1D,
    geometry::Geometry{Dim1},
    resolution,
    method::DiscretizationMethod=SimpleGrid(),
)
    ρ = discretize(geometry, resolution, :ρ, method)
    C11 = discretize(geometry, resolution, :C11, method)
    (ρ=ρ, C11=C11)
end

# 3D Full Vector EM (H-field formulation)
# LHS = curl × ε⁻¹ × curl, RHS = μ
function prepare_materials(
    ::FullVectorEM,
    geometry::Geometry{Dim3},
    resolution,
    method::DiscretizationMethod=SimpleGrid(),
)
    ε_inv = discretize(geometry, resolution, :ε_inv, method)
    μ = discretize(geometry, resolution, :μ, method)
    (ε_inv=ε_inv, μ=μ)
end

# 3D Full Elastic
# LHS = -∇·(C:∇), RHS = ρ
function prepare_materials(
    ::FullElastic,
    geometry::Geometry{Dim3},
    resolution,
    method::DiscretizationMethod=SimpleGrid(),
)
    ρ = discretize(geometry, resolution, :ρ, method)
    C11 = discretize(geometry, resolution, :C11, method)
    C12 = discretize(geometry, resolution, :C12, method)
    C44 = discretize(geometry, resolution, :C44, method)
    (ρ=ρ, C11=C11, C12=C12, C44=C44)
end

# ============================================================================
# Matrix construction
# ============================================================================

"""
    build_matrices(solver::Solver, k)

Build the eigenvalue problem matrices LHS * ψ = ω² * RHS * ψ for wave vector k.
"""
function build_matrices(solver::Solver{Dim2,TEWave}, k::AbstractVector{<:Real})
    basis = solver.basis
    mats = solver.material_arrays

    # Wave vector operators: K_x = k_x + G_x, K_y = k_y + G_y
    Kx = Diagonal([k[1] + G[1] for G in basis.G])
    Ky = Diagonal([k[2] + G[2] for G in basis.G])

    # LHS = Kx * ε⁻¹ * Kx + Ky * ε⁻¹ * Ky
    ε_inv_c = convolution_matrix(mats.ε_inv, basis)
    LHS = Kx * ε_inv_c * Kx + Ky * ε_inv_c * Ky

    # RHS = W (weight matrix)
    RHS = get_weight_matrix(solver)

    return LHS, RHS
end

function build_matrices(solver::Solver{Dim2,TMWave}, k::AbstractVector{<:Real})
    basis = solver.basis
    mats = solver.material_arrays

    Kx = Diagonal([k[1] + G[1] for G in basis.G])
    Ky = Diagonal([k[2] + G[2] for G in basis.G])

    # LHS = Kx * μ⁻¹ * Kx + Ky * μ⁻¹ * Ky
    μ_inv_c = convolution_matrix(mats.μ_inv, basis)
    LHS = Kx * μ_inv_c * Kx + Ky * μ_inv_c * Ky

    # RHS = W (weight matrix)
    RHS = get_weight_matrix(solver)

    return LHS, RHS
end

function build_matrices(solver::Solver{Dim2,SHWave}, k::AbstractVector{<:Real})
    basis = solver.basis
    mats = solver.material_arrays

    Kx = Diagonal([k[1] + G[1] for G in basis.G])
    Ky = Diagonal([k[2] + G[2] for G in basis.G])

    # LHS = Kx * C44 * Kx + Ky * C44 * Ky
    C44_c = convolution_matrix(mats.C44, basis)
    LHS = Kx * C44_c * Kx + Ky * C44_c * Ky

    # RHS = W (weight matrix)
    RHS = get_weight_matrix(solver)

    return LHS, RHS
end

function build_matrices(solver::Solver{Dim2,PSVWave}, k::AbstractVector{<:Real})
    basis = solver.basis
    mats = solver.material_arrays

    Kx = Diagonal([k[1] + G[1] for G in basis.G])
    Ky = Diagonal([k[2] + G[2] for G in basis.G])

    C11_c = convolution_matrix(mats.C11, basis)
    C12_c = convolution_matrix(mats.C12, basis)
    C44_c = convolution_matrix(mats.C44, basis)

    # Block matrix construction
    # K_xx = ∂_x C₁₁ ∂_x + ∂_y C₄₄ ∂_y
    K_xx = Kx * C11_c * Kx + Ky * C44_c * Ky

    # K_yy = ∂_x C₄₄ ∂_x + ∂_y C₁₁ ∂_y
    K_yy = Kx * C44_c * Kx + Ky * C11_c * Ky

    # K_xy = ∂_x C₁₂ ∂_y + ∂_y C₄₄ ∂_x
    K_xy = Kx * C12_c * Ky + Ky * C44_c * Kx

    # K_yx = ∂_x C₄₄ ∂_y + ∂_y C₁₂ ∂_x
    K_yx = Kx * C44_c * Ky + Ky * C12_c * Kx

    # Assemble 2N × 2N block matrices
    LHS = [K_xx K_xy; K_yx K_yy]

    # RHS = W (weight matrix)
    RHS = get_weight_matrix(solver)

    return LHS, RHS
end

# ============================================================================
# 1D Matrix construction
# ============================================================================

function build_matrices(solver::Solver{Dim1,Photonic1D}, k::Real)
    basis = solver.basis
    mats = solver.material_arrays

    # Wave vector operator: K = k + G
    K = Diagonal([k + G[1] for G in basis.G])

    # LHS = K * ε⁻¹ * K
    ε_inv_c = convolution_matrix(mats.ε_inv, basis)
    LHS = K * ε_inv_c * K

    # RHS = W (weight matrix)
    RHS = get_weight_matrix(solver)

    return LHS, RHS
end

function build_matrices(solver::Solver{Dim1,Longitudinal1D}, k::Real)
    basis = solver.basis
    mats = solver.material_arrays

    # Wave vector operator: K = k + G
    K = Diagonal([k + G[1] for G in basis.G])

    # LHS = K * C11 * K
    C11_c = convolution_matrix(mats.C11, basis)
    LHS = K * C11_c * K

    # RHS = W (weight matrix)
    RHS = get_weight_matrix(solver)

    return LHS, RHS
end

# 1D convenience methods accepting Vector{Float64} (for compatibility with solve interface)
function build_matrices(solver::Solver{Dim1,Photonic1D}, k::Vector{Float64})
    build_matrices(solver, k[1])
end
function build_matrices(solver::Solver{Dim1,Longitudinal1D}, k::Vector{Float64})
    build_matrices(solver, k[1])
end

# ============================================================================
# 3D Matrix construction
# ============================================================================

"""
    skew_matrix(k)

Create 3x3 skew-symmetric matrix for cross product: [k]× such that [k]× v = k × v
"""
function skew_matrix(k::AbstractVector)
    [
        0.0 -k[3] k[2];
        k[3] 0.0 -k[1];
        -k[2] k[1] 0.0
    ]
end

"""
    build_matrices(solver::Solver{Dim3, FullVectorEM}, k)

Build 3N×3N matrices for 3D photonic crystal (H-field formulation).

Formulation: curl × ε⁻¹ × curl H = (ω²/c²) μ H

Storage order: [H_x; H_y; H_z] (component-order for FFT efficiency)
"""
function build_matrices(solver::Solver{Dim3,FullVectorEM}, k::AbstractVector{<:Real})
    basis = solver.basis
    mats = solver.material_arrays

    # Convolution matrix for LHS
    ε_inv_c = convolution_matrix(mats.ε_inv, basis)

    # Wave vectors K = k + G for each plane wave
    Kx = Diagonal([k[1] + G[1] for G in basis.G])
    Ky = Diagonal([k[2] + G[2] for G in basis.G])
    Kz = Diagonal([k[3] + G[3] for G in basis.G])

    # Build LHS: 9 blocks L_ij where L_ij = curl_i × ε⁻¹ × curl_j
    # curl_i is the i-th row of the curl operator in Fourier space
    # curl = i(k+G)× = i * skew_matrix(K)
    # For component i of curl×(ε⁻¹×(curl×H)):
    #   (curl × (ε⁻¹ × (curl × H)))_x = Ky*(ε⁻¹*(Kx*Hy - Ky*Hx)) - Kz*(ε⁻¹*(Kz*Hx - Kx*Hz))
    #                                 = -Ky*ε⁻¹*Ky*Hx + Ky*ε⁻¹*Kx*Hy - Kz*ε⁻¹*Kz*Hx + Kz*ε⁻¹*Kx*Hz
    #                                 = -(Ky*ε⁻¹*Ky + Kz*ε⁻¹*Kz)*Hx + Ky*ε⁻¹*Kx*Hy + Kz*ε⁻¹*Kx*Hz

    # L_xx = -(Ky² + Kz²)*ε⁻¹ (but with convolution in between)
    # More precisely: L_xx = Ky*ε⁻¹*Ky + Kz*ε⁻¹*Kz
    L_xx = Ky * ε_inv_c * Ky + Kz * ε_inv_c * Kz
    L_yy = Kx * ε_inv_c * Kx + Kz * ε_inv_c * Kz
    L_zz = Kx * ε_inv_c * Kx + Ky * ε_inv_c * Ky

    # Off-diagonal blocks
    L_xy = -Ky * ε_inv_c * Kx
    L_xz = -Kz * ε_inv_c * Kx
    L_yx = -Kx * ε_inv_c * Ky
    L_yz = -Kz * ε_inv_c * Ky
    L_zx = -Kx * ε_inv_c * Kz
    L_zy = -Ky * ε_inv_c * Kz

    # Assemble 3N × 3N block matrix
    LHS = [
        L_xx L_xy L_xz;
        L_yx L_yy L_yz;
        L_zx L_zy L_zz
    ]

    # RHS = W (weight matrix)
    RHS = get_weight_matrix(solver)

    return LHS, RHS
end

"""
    build_matrices(solver::Solver{Dim3, FullElastic}, k)

Build 3N×3N matrices for 3D phononic crystal.

Formulation: -∇·(C:∇u) = ω² ρ u

Storage order: [u_x; u_y; u_z] (component-order)
"""
function build_matrices(solver::Solver{Dim3,FullElastic}, k::AbstractVector{<:Real})
    basis = solver.basis
    mats = solver.material_arrays

    # Convolution matrices for LHS
    C11_c = convolution_matrix(mats.C11, basis)
    C12_c = convolution_matrix(mats.C12, basis)
    C44_c = convolution_matrix(mats.C44, basis)

    # Wave vectors
    Kx = Diagonal([k[1] + G[1] for G in basis.G])
    Ky = Diagonal([k[2] + G[2] for G in basis.G])
    Kz = Diagonal([k[3] + G[3] for G in basis.G])

    # Build stiffness matrix blocks for isotropic material
    # K_xx = Kx*C11*Kx + Ky*C44*Ky + Kz*C44*Kz
    K_xx = Kx * C11_c * Kx + Ky * C44_c * Ky + Kz * C44_c * Kz
    K_yy = Kx * C44_c * Kx + Ky * C11_c * Ky + Kz * C44_c * Kz
    K_zz = Kx * C44_c * Kx + Ky * C44_c * Ky + Kz * C11_c * Kz

    # Off-diagonal blocks
    K_xy = Kx * C12_c * Ky + Ky * C44_c * Kx
    K_xz = Kx * C12_c * Kz + Kz * C44_c * Kx
    K_yx = Ky * C12_c * Kx + Kx * C44_c * Ky
    K_yz = Ky * C12_c * Kz + Kz * C44_c * Ky
    K_zx = Kz * C12_c * Kx + Kx * C44_c * Kz
    K_zy = Kz * C12_c * Ky + Ky * C44_c * Kz

    # Assemble 3N × 3N block matrix
    LHS = [
        K_xx K_xy K_xz;
        K_yx K_yy K_yz;
        K_zx K_zy K_zz
    ]

    # RHS = W (weight matrix)
    RHS = get_weight_matrix(solver)

    return LHS, RHS
end

# ============================================================================
# get_weight_matrix: Return weight matrix for inner products
# ============================================================================

"""
    get_weight_matrix(solver::Solver)

Return the weight matrix W for computing inner products of eigenvectors.

The weight matrix is the RHS matrix of the generalized eigenvalue problem
(LHS * v = ω² * RHS * v). It defines the inner product in the function space:
for eigenvectors v₁ and v₂, the overlap is computed as v₁' * W * v₂.

For normalized eigenvectors from `solve_at_k_with_vectors`, we have:
`vecs' * W * vecs ≈ I`

# Weight matrices by wave type:
| WaveType | W |
|----------|---|
| TE | μ (permeability) |
| TM | ε (permittivity) |
| SH | ρ (density) |
| PSV | diag(ρ, ρ) |
| FullVectorEM | diag(ε, ε, ε) |
| FullElastic | diag(ρ, ρ, ρ) |
| Photonic1D | ε |
| Longitudinal1D | ρ |

# Returns
- `W::Matrix{ComplexF64}`: Weight matrix

# Example
```julia
W = get_weight_matrix(solver)
ω, vecs = solve_at_k_with_vectors(solver, k, DenseMethod(); bands=1:4)
overlap = vecs' * W * vecs  # ≈ I(4)
```

See also: [`solve_at_k_with_vectors`](@ref), [`build_matrices`](@ref)
"""
function get_weight_matrix(solver::Solver{Dim2,TEWave})
    # TE: W = μ (permeability convolution matrix)
    return convolution_matrix(solver.material_arrays.μ, solver.basis)
end

function get_weight_matrix(solver::Solver{Dim2,TMWave})
    # TM: W = ε (permittivity convolution matrix)
    return convolution_matrix(solver.material_arrays.ε, solver.basis)
end

function get_weight_matrix(solver::Solver{Dim2,SHWave})
    # SH: W = ρ (density convolution matrix)
    return convolution_matrix(solver.material_arrays.ρ, solver.basis)
end

function get_weight_matrix(solver::Solver{Dim2,PSVWave})
    # PSV: W = diag(ρ, ρ)
    N = solver.basis.num_pw
    ρ_c = convolution_matrix(solver.material_arrays.ρ, solver.basis)
    Z = zeros(ComplexF64, N, N)
    return [
        ρ_c Z
        Z ρ_c
    ]
end

function get_weight_matrix(solver::Solver{Dim1,Photonic1D})
    # 1D Photonic: W = μ (consistent with RHS in build_matrices)
    return convolution_matrix(solver.material_arrays.μ, solver.basis)
end

function get_weight_matrix(solver::Solver{Dim1,Longitudinal1D})
    # 1D Elastic: W = ρ
    return convolution_matrix(solver.material_arrays.ρ, solver.basis)
end

function get_weight_matrix(solver::Solver{Dim3,FullVectorEM})
    # 3D EM: W = diag(μ, μ, μ) (consistent with RHS in build_matrices)
    N = solver.basis.num_pw
    μ_c = convolution_matrix(solver.material_arrays.μ, solver.basis)
    Z = zeros(ComplexF64, N, N)
    return [
        μ_c Z Z
        Z μ_c Z
        Z Z μ_c
    ]
end

function get_weight_matrix(solver::Solver{Dim3,FullElastic})
    # 3D Elastic: W = diag(ρ, ρ, ρ)
    N = solver.basis.num_pw
    ρ_c = convolution_matrix(solver.material_arrays.ρ, solver.basis)
    Z = zeros(ComplexF64, N, N)
    return [
        ρ_c Z Z
        Z ρ_c Z
        Z Z ρ_c
    ]
end

# ============================================================================
# Solve eigenvalue problem
# ============================================================================

"""
    solve(solver::Solver, k; bands=1:10)

Solve the eigenvalue problem at wave vector k.

# Arguments
- `solver`: The solver
- `k`: Wave vector (in units of reciprocal lattice vectors or absolute)
- `bands`: Which bands to return (default: 1:10)

# Returns
- `frequencies`: Eigenfrequencies (sorted, positive)
- `modes`: Corresponding eigenvectors
"""
function solve(solver::Solver, k; bands=1:10)
    solve_impl(solver, k, solver.method; bands=bands)
end

# Convenience method with SVector
function solve(solver::Solver{Dim2}, k::SVector{2}; bands=1:10)
    solve(solver, Vector(k); bands=bands)
end

# Convenience method with tuple
function solve(solver::Solver{Dim2}, k::Tuple{Real,Real}; bands=1:10)
    solve(solver, [Float64(k[1]), Float64(k[2])]; bands=bands)
end

# Convenience methods for 3D
function solve(solver::Solver{Dim3}, k::SVector{3}; bands=1:10)
    solve(solver, Vector(k); bands=bands)
end

function solve(solver::Solver{Dim3}, k::Tuple{Real,Real,Real}; bands=1:10)
    solve(solver, [Float64(k[1]), Float64(k[2]), Float64(k[3])]; bands=bands)
end

# ============================================================================
# solve_at_k: Single k-point solver with manual control
# ============================================================================

"""
    solve_at_k(solver, k, method; bands=1:10, X0=nothing, P=nothing)

Solve eigenvalue problem at a single k-point with explicit control over method,
initial vectors, and preconditioner. Returns only frequencies.

For eigenvectors, use [`solve_at_k_with_vectors`](@ref) instead.

# Arguments
- `solver::Solver`: Solver instance
- `k`: Wave vector (2D: Vector{Float64}, 1D: Float64)
- `method::SolverMethod`: Solver method (DenseMethod, LOBPCGMethod, etc.)

# Keyword Arguments
- `bands`: Range of bands to compute (default: 1:10)
- `X0`: Initial guess for eigenvectors (default: nothing = auto)
  - Size: `dim × nev` where `dim = matrix_dimension(solver)` and `nev = maximum(bands)`
  - Type: `Matrix{ComplexF64}`
  - If `nothing`, random orthonormal vectors are used
- `P`: Preconditioner (default: nothing = use method's preconditioner setting)
  - Must implement `ldiv!(y, P, x)`

# Returns
- `Vector{Float64}`: Frequencies for the requested bands

# Example
```julia
solver = Solver(PSVWave(), geo, (64, 64); cutoff=20)
dim = matrix_dimension(solver)  # 2514 for cutoff=20

# Frequencies only
freqs = solve_at_k(solver, k, LOBPCGMethod(); bands=1:20)

# With eigenvectors - use solve_at_k_with_vectors
freqs, vecs = solve_at_k_with_vectors(solver, k, DenseMethod(); bands=1:20)

# Custom preconditioner
using LinearAlgebra
LHS, _ = build_matrices(solver, k)
P = Diagonal(1.0 ./ diag(LHS))
freqs = solve_at_k(solver, k, LOBPCGMethod(); bands=1:20, P=P)
```

See also: [`solve_at_k_with_vectors`](@ref), [`solve`](@ref), [`matrix_dimension`](@ref)
"""
function solve_at_k(
    solver::Solver, k, method::SolverMethod; bands=1:10, X0=nothing, P=nothing
)
    frequencies, _ = _solve_at_k_impl(solver, k, method; bands=bands, X0=X0, P=P)
    return frequencies
end

"""
    solve_at_k_with_vectors(solver, k, method; bands=1:10, X0=nothing, P=nothing)

Solve eigenvalue problem at a single k-point and return both frequencies and eigenvectors.

This function always returns eigenvectors, unlike `solve_at_k` which returns only frequencies
by default. Use this when you need the eigenvectors for further analysis (e.g., computing
overlaps, mode profiles, or topological invariants).

# Arguments
- `solver::AbstractSolver`: The solver object
- `k`: Wave vector (Real for 1D, AbstractVector for 2D/3D)
- `method::SolverMethod`: Solver method (DenseMethod(), KrylovKitMethod(), LOBPCGMethod())

# Keyword Arguments
- `bands`: Range of bands to compute (default: 1:10)
- `X0`: Initial guess for eigenvectors (for LOBPCG, optional)
- `P`: Preconditioner (for LOBPCG, optional)

# Returns
- `(frequencies, eigenvectors)`: Tuple of frequencies (Vector{Float64}) and eigenvectors (Matrix{ComplexF64})

# Example
```julia
ω, vecs = solve_at_k_with_vectors(solver, [0.1, 0.2], DenseMethod(); bands=1:4)
W = get_weight_matrix(solver)
overlap = vecs' * W * vecs  # Should be ≈ I (orthonormal)
```

See also: [`solve_at_k`](@ref), [`get_weight_matrix`](@ref), [`build_matrices`](@ref)
"""
function solve_at_k_with_vectors(
    solver::Solver, k, method::SolverMethod; bands=1:10, X0=nothing, P=nothing
)
    frequencies, eigenvectors = _solve_at_k_impl(solver, k, method; bands=bands, X0=X0, P=P)
    return (frequencies, eigenvectors)
end

# Internal implementation dispatching by method type
function _solve_at_k_impl(
    solver::Solver, k, method::DenseMethod; bands=1:10, X0=nothing, P=nothing
)
    # Dense method ignores X0 and P
    LHS, RHS = build_matrices(solver, k)
    _solve_dense(LHS, RHS, bands, method.shift)
end

function _solve_at_k_impl(
    solver::Solver, k, method::KrylovKitMethod; bands=1:10, X0=nothing, P=nothing
)
    # KrylovKit currently ignores X0 and P
    _solve_krylovkit(solver, k, method, bands)
end

function _solve_at_k_impl(
    solver::Solver, k, method::LOBPCGMethod; bands=1:10, X0=nothing, P=nothing
)
    _solve_lobpcg_at_k(solver, k, method; bands=bands, X0=X0, P=P)
end

# Fallback
function _solve_at_k_impl(
    solver::Solver, k, method::SolverMethod; bands=1:10, X0=nothing, P=nothing
)
    error("solve_at_k not implemented for method: $(typeof(method))")
end

# ============================================================================
# Dense method implementation (2D)
# ============================================================================

function solve_impl(
    solver::Solver{Dim2}, k::Vector{<:Real}, method::DenseMethod; bands=1:10
)
    LHS, RHS = build_matrices(solver, k)
    _solve_dense(LHS, RHS, bands, method.shift)
end

# ============================================================================
# Dense method implementation (3D)
# ============================================================================

function solve_impl(
    solver::Solver{Dim3}, k::Vector{<:Real}, method::DenseMethod; bands=1:10
)
    LHS, RHS = build_matrices(solver, k)
    _solve_dense(LHS, RHS, bands, method.shift)
end

# ============================================================================
# Dense method implementation (1D)
# ============================================================================

function solve_impl(solver::Solver{Dim1}, k::Real, method::DenseMethod; bands=1:10)
    LHS, RHS = build_matrices(solver, k)
    _solve_dense(LHS, RHS, bands, method.shift)
end

# ============================================================================
# Dense solver core (type-stable inner function)
# ============================================================================

"""
    _solve_dense(LHS, RHS, bands, shift=0.0)

Core dense eigenvalue solver. Type-stable implementation.

# Arguments
- `LHS`: Left-hand side matrix (stiffness-like)
- `RHS`: Right-hand side matrix (mass-like)
- `bands`: Range of bands to return
- `shift`: Minimum eigenvalue cutoff (eigenvalues with ω² < shift are filtered out)
"""
function _solve_dense(LHS::AbstractMatrix, RHS::AbstractMatrix, bands, shift::Real=0.0)
    # Solve generalized eigenvalue problem
    # Handle potential singularity at Γ point
    eigenvalues, eigenvectors = try
        eigen(Hermitian(Array(LHS)), Hermitian(Array(RHS)))
    catch
        # Fallback: standard eigenvalue problem
        vals, vecs = eigen(Array(RHS) \ Array(LHS))
        vals, vecs
    end

    # Convert to frequencies (ω = sqrt(eigenvalue))
    # For photonic: eigenvalue = (ω/c)², for phononic: eigenvalue = ω²
    ω² = real.(eigenvalues)

    # Handle small negative eigenvalues (numerical noise)
    ω² = max.(ω², 0.0)
    frequencies = sqrt.(ω²)

    # Sort by frequency
    perm = sortperm(frequencies)
    frequencies = frequencies[perm]
    eigenvectors = eigenvectors[:, perm]

    # Filter out eigenvalues below shift (spurious modes at ω ≈ 0)
    if shift > 0
        valid_mask = (frequencies .^ 2) .>= shift
        valid_idx = findall(valid_mask)
        if !isempty(valid_idx)
            frequencies = frequencies[valid_idx]
            eigenvectors = eigenvectors[:, valid_idx]
        end
    end

    # Select requested bands
    return _select_bands(frequencies, eigenvectors, bands)
end

# ============================================================================
# KrylovKit iterative method implementation (2D)
# ============================================================================

function solve_impl(
    solver::Solver{Dim2,W}, k::Vector{<:Real}, method::KrylovKitMethod; bands=1:10
) where {W<:WaveType}
    _solve_krylovkit(solver, k, method, bands)
end

# ============================================================================
# KrylovKit iterative method implementation (1D)
# ============================================================================

function solve_impl(
    solver::Solver{Dim1,W}, k::Real, method::KrylovKitMethod; bands=1:10
) where {W<:WaveType}
    _solve_krylovkit(solver, [Float64(k)], method, bands)
end

# ============================================================================
# KrylovKit iterative method implementation (3D)
# ============================================================================

function solve_impl(
    solver::Solver{Dim3,W}, k::Vector{<:Real}, method::KrylovKitMethod; bands=1:10
) where {W<:WaveType}
    _solve_krylovkit(solver, Vector{Float64}(k), method, bands)
end

# ============================================================================
# KrylovKit solver core
# ============================================================================

"""
    _solve_krylovkit(solver, k, method, bands)

Core iterative eigenvalue solver using KrylovKit.jl.

When `method.shift > 0`, uses shift-and-invert spectral transformation:
- Original: `A x = λ B x`
- Transformed: `(A - σB)⁻¹ B x = μ x` where `μ = 1/(λ - σ)`

This finds eigenvalues closest to σ, which is useful for:
- Skipping spurious longitudinal modes in 3D H-field formulation
- Targeting specific frequency ranges
"""
function _solve_krylovkit(
    solver::Solver{D,W}, k::Vector{Float64}, method::KrylovKitMethod, bands
) where {D,W}
    σ = method.shift

    if σ > 0
        if method.matrix_free
            # Shift-and-invert mode: matrix-free with iterative inner solver
            _solve_krylovkit_shifted_matrixfree(solver, k, method, bands)
        else
            # Shift-and-invert mode: use dense matrices with factorization
            _solve_krylovkit_shifted(solver, k, method, bands)
        end
    else
        # Standard mode: use matrix-free operators
        _solve_krylovkit_standard(solver, k, method, bands)
    end
end

"""
    _solve_krylovkit_standard(solver, k, method, bands)

Standard generalized eigenvalue solver (no shift).
Uses matrix-free operators for memory efficiency.

For phononic problems, applies scaling to normalize eigenvalues to O(1),
which improves numerical stability in KrylovKit.geneigsolve.
"""
function _solve_krylovkit_standard(
    solver::Solver{D,W}, k::Vector{Float64}, method::KrylovKitMethod, bands
) where {D,W}
    # Create matrix-free operator
    op = MatrixFreeOperator(solver, k)
    N = solver.basis.num_pw
    nc = ncomponents(solver.wave)
    dim = N * nc

    # Number of eigenvalues to compute
    nev = _nev(bands, dim)

    # Estimate scaling factor for phononic problems
    # For elastic waves, eigenvalues ω² ~ (stiffness/density) * k²
    # We scale so that eigenvalues become O(1)
    scale_factor = _estimate_eigenvalue_scale(solver, k)

    # Define functions for scaled A*x and B*x
    # Original: A x = λ B x
    # Scaled: (A/s) x = (λ/s) B x, where s = scale_factor
    function Ax(x)
        y = zeros(ComplexF64, dim)
        apply_lhs!(y, op, x)
        return y ./ scale_factor
    end

    function Bx(x)
        y = zeros(ComplexF64, dim)
        apply_rhs!(y, op, x)
        return y
    end

    # Initial vector
    x0 = randn(ComplexF64, dim)
    x0 ./= norm(x0)

    # Solve generalized eigenvalue problem: (A/s)*x = λ̃*B*x
    # where λ̃ = λ/s (scaled eigenvalue)
    vals, vecs, info = KrylovKit.geneigsolve(
        (Ax, Bx),
        x0,
        nev,
        :SR;
        tol=method.tol,
        maxiter=method.maxiter,
        krylovdim=method.krylovdim,
        verbosity=method.verbosity,
        issymmetric=true,
        isposdef=true,
    )

    # Check convergence
    if info.converged < nev && method.verbosity > 0
        @warn "KrylovKit: Only $(info.converged) of $nev eigenvalues converged"
    end

    # Convert back to original eigenvalues: λ = λ̃ * s
    ω² = real.(vals) .* scale_factor
    ω² = max.(ω², 0.0)  # Handle small negative values (numerical noise)
    frequencies = sqrt.(ω²)

    # Build eigenvector matrix
    eigenvectors = stack(vecs)

    # Sort by frequency (should already be sorted, but ensure)
    perm = sortperm(frequencies)
    frequencies = frequencies[perm]
    eigenvectors = eigenvectors[:, perm]

    # Select requested bands
    return _select_bands(frequencies, eigenvectors, bands)
end

"""
    _estimate_eigenvalue_scale(solver, k)

Estimate a characteristic scale for eigenvalues to normalize them to O(1).
This improves numerical stability in iterative eigensolvers.

For photonic: ω² ~ (c/a)² ~ 1 (already normalized)
For phononic: ω² ~ (c_s * k)² where c_s is shear wave speed
"""
function _estimate_eigenvalue_scale(solver::Solver{D,W}, k::Vector{Float64}) where {D,W}
    _estimate_eigenvalue_scale(solver.wave, solver, k)
end

# Photonic: eigenvalues are already O(1) due to normalization
_estimate_eigenvalue_scale(::PhotonicWave, solver, k) = 1.0

# Phononic: estimate from material properties and k (dispatch on material type)
function _estimate_eigenvalue_scale(::PhononicWave, solver, k)
    _eigenvalue_scale_from_material(solver.geometry.background, k)
end

# Multiple dispatch for material types
function _eigenvalue_scale_from_material(mat::IsotropicElastic, k)
    # Characteristic wave speed squared: c² = C11/ρ = (λ + 2μ) / ρ for P-wave
    c² = mat.C11 / mat.ρ
    k_mag = norm(k)
    # Eigenvalue scale: ω² ~ c² * k²
    max(c² * k_mag^2, 1.0)  # Ensure at least 1
end

# Fallback for other material types
_eigenvalue_scale_from_material(::Material, k) = 1.0

"""
    _solve_krylovkit_shifted(solver, k, method, bands)

Shift-and-invert eigenvalue solver.
Transforms `A x = λ B x` to `(A - σB)⁻¹ B x = μ x` where `μ = 1/(λ - σ)`.

Finds eigenvalues closest to σ from above (λ > σ), useful for:
- Skipping spurious modes (e.g., longitudinal modes in 3D with λ ≈ 0)
- Targeting specific frequency ranges

# Algorithm
Uses KrylovKit.eigsolve with `:LR` (largest real) to find eigenvalues μ with
the largest real part. Since μ = 1/(λ - σ):
- λ > σ ⟹ μ > 0 (positive, selected by :LR)
- λ < σ ⟹ μ < 0 (negative, filtered out)

The positive μ values are explicitly filtered to ensure only eigenvalues
above the shift are returned. This provides robustness when:
- The user specifies an inappropriate σ (larger than target eigenvalues)
- KrylovKit returns more eigenvalues than requested

# Note
The shift σ acts as a lower bound: only eigenvalues λ > σ are returned.
For 3D photonic crystals, σ = 0.01 effectively skips longitudinal modes (λ ≈ 0).
"""
function _solve_krylovkit_shifted(
    solver::Solver{D,W}, k::Vector{Float64}, method::KrylovKitMethod, bands
) where {D,W}
    σ = method.shift

    # Build dense matrices (required for factorization)
    LHS, RHS = build_matrices(solver, k)
    dim = size(LHS, 1)

    # Number of eigenvalues to compute
    nev = _nev(bands, dim)

    # Compute shifted matrix and factorize: (A - σB)
    shifted = LHS - σ * RHS
    shifted_fact = lu(shifted)

    # Define operator: (A - σB)⁻¹ B x
    function shifted_op(x)
        shifted_fact \ (RHS * x)
    end

    # Initial vector
    x0 = randn(ComplexF64, dim)
    x0 ./= norm(x0)

    # For eigenvalues λ > σ, we have μ = 1/(λ-σ) > 0
    # For eigenvalues λ < σ (including λ ≈ 0), we have μ < 0
    # We want positive μ (largest real), which gives λ closest to σ from above
    # Request more eigenvalues to ensure we get enough positive ones
    nev_request = min(dim, nev * 3 + 10)

    μ_vals, vecs, info = KrylovKit.eigsolve(
        shifted_op,
        x0,
        nev_request,
        :LR;  # Largest Real (positive μ)
        tol=method.tol,
        maxiter=method.maxiter,
        krylovdim=max(method.krylovdim, nev_request + 10),
        verbosity=method.verbosity,
        ishermitian=true,
    )

    # Filter for positive μ (corresponding to λ > σ)
    positive_indices = findall(x -> real(x) > 0, μ_vals)

    if isempty(positive_indices)
        @warn "No eigenvalues found above shift σ = $σ. Try a smaller shift."
        # Fall back to returning what we have
        positive_indices = 1:min(nev, length(μ_vals))
    end

    μ_positive = μ_vals[positive_indices]
    vecs_positive = vecs[positive_indices]

    # Check convergence
    if length(positive_indices) < nev && method.verbosity > 0
        @warn "KrylovKit (shifted): Only $(length(positive_indices)) eigenvalues found above shift"
    end

    # Convert back to original eigenvalues: λ = σ + 1/μ
    λ_vals = σ .+ 1.0 ./ real.(μ_positive)

    # Convert to frequencies
    ω² = max.(λ_vals, 0.0)  # Handle small negative values
    frequencies = sqrt.(ω²)

    # Build eigenvector matrix
    eigenvectors = stack(vecs_positive)

    # Sort by frequency (smallest first)
    perm = sortperm(frequencies)
    frequencies = frequencies[perm]
    eigenvectors = eigenvectors[:, perm]

    # Select requested bands (uses multiple dispatch)
    return _select_first_bands(frequencies, eigenvectors, bands)
end

"""
    _solve_krylovkit_shifted_matrixfree(solver, k, method, bands)

Matrix-free shift-and-invert eigenvalue solver.

Unlike `_solve_krylovkit_shifted`, this version uses O(N) memory instead of O(N²)
by avoiding explicit construction of dense matrices. Instead, it uses:
1. Matrix-free operators (apply_lhs!, apply_rhs!) via FFT
2. Iterative inner solver (Krylov.cg) to solve (A - σB)y = Bx

This is particularly useful for large 3D calculations where dense matrix storage
becomes prohibitive.

# Algorithm
1. Create matrix-free ShiftedOperator S = (A - σB)
2. For each application of (A - σB)⁻¹ B x:
   a. Compute z = B x using apply_rhs!
   b. Solve S y = z using iterative CG
3. Use KrylovKit.eigsolve with :LR to find largest μ = 1/(λ - σ)
4. Convert back: λ = σ + 1/μ

# Memory Complexity
| Component | Dense | Matrix-Free |
|-----------|-------|-------------|
| LHS, RHS  | O(N²) | O(N)        |
| Shifted   | O(N²) | O(N)        |
| LU factor | O(N²) | N/A         |
| Total     | ~3N²  | ~4N         |
"""
function _solve_krylovkit_shifted_matrixfree(
    solver::Solver{D,W}, k::Vector{Float64}, method::KrylovKitMethod, bands
) where {D,W}
    σ = method.shift

    # Create matrix-free operator
    op = MatrixFreeOperator(solver, k)
    N = solver.basis.num_pw
    nc = ncomponents(solver.wave)
    dim = N * nc

    # Number of eigenvalues to compute
    nev = _nev(bands, dim)

    # Create shifted operator S = (A - σB)
    S = ShiftedOperator(op, σ)

    # Temporary vectors for B*x computation
    z = zeros(ComplexF64, dim)

    # Inner solver tolerance (tighter than outer)
    inner_tol = method.tol * 0.1
    inner_maxiter = 500

    # Function: (A - σB)⁻¹ B x
    function shifted_inv_B(x)
        # z = B * x
        apply_rhs!(z, op, x)

        # Solve (A - σB) y = z using CG
        # Note: Krylov.cg expects real symmetric or complex Hermitian
        y, stats = Krylov.cg(S, z; atol=inner_tol, rtol=inner_tol, itmax=inner_maxiter)

        if !stats.solved && method.verbosity > 0
            @warn "Inner CG solver did not converge: $(stats.status)"
        end

        return y
    end

    # Initial vector
    x0 = randn(ComplexF64, dim)
    x0 ./= norm(x0)

    # Request more eigenvalues to ensure we get enough positive ones
    nev_request = min(dim, nev * 3 + 10)

    # For eigenvalues λ > σ, we have μ = 1/(λ-σ) > 0
    # For eigenvalues λ < σ, we have μ < 0
    # We want positive μ (largest real), which gives λ closest to σ from above
    μ_vals, vecs, info = KrylovKit.eigsolve(
        shifted_inv_B,
        x0,
        nev_request,
        :LR;  # Largest Real (positive μ)
        tol=method.tol,
        maxiter=method.maxiter,
        krylovdim=max(method.krylovdim, nev_request + 10),
        verbosity=method.verbosity,
        ishermitian=true,
    )

    # Filter for positive μ (corresponding to λ > σ)
    positive_indices = findall(x -> real(x) > 0, μ_vals)

    if isempty(positive_indices)
        @warn "Matrix-free shifted: No eigenvalues found above shift σ = $σ. Try a smaller shift."
        # Fall back to returning what we have
        positive_indices = 1:min(nev, length(μ_vals))
    end

    μ_positive = μ_vals[positive_indices]
    vecs_positive = vecs[positive_indices]

    # Check convergence
    if length(positive_indices) < nev && method.verbosity > 0
        @warn "Matrix-free shifted: Only $(length(positive_indices)) eigenvalues found above shift"
    end

    # Convert back to original eigenvalues: λ = σ + 1/μ
    λ_vals = σ .+ 1.0 ./ real.(μ_positive)

    # Convert to frequencies
    ω² = max.(λ_vals, 0.0)  # Handle small negative values
    frequencies = sqrt.(ω²)

    # Build eigenvector matrix
    eigenvectors = stack(vecs_positive)

    # Sort by frequency (smallest first)
    perm = sortperm(frequencies)
    frequencies = frequencies[perm]
    eigenvectors = eigenvectors[:, perm]

    # Select requested bands (uses multiple dispatch)
    return _select_first_bands(frequencies, eigenvectors, bands)
end

# ============================================================================
# LOBPCG method implementation
# ============================================================================

"""
    solve_impl(solver::Solver, k, method::LOBPCGMethod; bands=1:10)

Solve eigenvalue problem using LOBPCG (Locally Optimal Block Preconditioned
Conjugate Gradient) method from IterativeSolvers.jl.

Based on Knyazev (2001), SIAM J. Sci. Comput. Vol.23, No.2, pp.517-541.

LOBPCG is particularly effective for:
- Symmetric generalized eigenvalue problems A x = λ B x
- Problems where B is positive definite (mass matrix)
- Computing multiple smallest eigenvalues simultaneously (block method)

Unlike KrylovKit, LOBPCG does not require explicit eigenvalue scaling for
phononic problems, as it works directly with the symmetric structure.

When `method.shift > 0`, uses shift-and-invert spectral transformation:
- Original: `A x = λ B x`
- Transformed: `(A - σB)⁻¹ B x = μ x` where `μ = 1/(λ - σ)`
- LOBPCG finds largest μ, corresponding to smallest λ > σ
"""
function _solve_lobpcg(solver::Solver, k, method::LOBPCGMethod; bands=1:10)
    σ = method.shift

    if σ > 0
        _solve_lobpcg_shifted(solver, k, method; bands=bands)
    else
        _solve_lobpcg_standard(solver, k, method; bands=bands)
    end
end

"""
    _solve_lobpcg_at_k(solver, k, method; bands, X0, P)

LOBPCG solver with explicit initial vectors and preconditioner support.
Used by solve_at_k for fine-grained control.

Supports:
- X0: Initial eigenvector guess (warm start)
- P: Custom preconditioner (overrides method.preconditioner)
- scale: Matrix scaling for better conditioning
"""
function _solve_lobpcg_at_k(
    solver::Solver, k, method::LOBPCGMethod; bands=1:10, X0=nothing, P=nothing
)
    # Note: shift-and-invert not yet supported with X0/P
    if method.shift > 0
        return _solve_lobpcg_shifted(solver, k, method; bands=bands)
    end

    # Build dense matrices
    LHS, RHS = build_matrices(solver, k)
    dim = size(LHS, 1)
    nev = _nev(bands, dim)

    # Apply scaling if enabled
    scale_factor = 1.0
    if method.scale
        scale_factor = maximum(abs.(LHS))
        A = Hermitian(LHS / scale_factor)
    else
        A = Hermitian(LHS)
    end
    B = Hermitian(RHS)

    # Prepare initial vectors
    if X0 === nothing
        # Random orthonormal vectors
        X0_use = randn(ComplexF64, dim, nev)
        X0_use, _ = qr(X0_use)
        X0_use = Matrix(X0_use)
    else
        # Use provided initial vectors (re-orthonormalize)
        X0_use, _ = qr(X0)
        X0_use = Matrix(X0_use)
    end

    # Prepare preconditioner
    if P !== nothing
        # Use provided preconditioner
        P_use = P
    elseif method.preconditioner == :diagonal
        # Diagonal preconditioner
        d = diag(A)
        d_safe = [abs(x) > 1e-10 ? x : 1e-10 for x in d]
        P_use = Diagonal(1.0 ./ d_safe)
    elseif method.preconditioner == :none
        P_use = nothing
    else
        # Assume it's a custom preconditioner object
        P_use = method.preconditioner
    end

    # Solve using LOBPCG
    results = IterativeSolvers.lobpcg(
        A, B, false, X0_use; P=P_use, tol=method.tol, maxiter=method.maxiter
    )

    # Extract eigenvalues and eigenvectors
    λ_vals = results.λ
    eigvecs = results.X

    # Rescale eigenvalues if scaling was applied
    if method.scale
        λ_vals = λ_vals * scale_factor
    end

    # Convert to frequencies: ω = √λ
    ω² = real.(λ_vals)
    ω² = max.(ω², 0.0)
    frequencies = sqrt.(ω²)

    # Sort by frequency
    perm = sortperm(frequencies)
    frequencies = frequencies[perm]
    eigenvectors = eigvecs[:, perm]

    # Select requested bands
    return _select_bands(frequencies, eigenvectors, bands)
end

"""
    _solve_lobpcg_standard(solver, k, method; bands)

Standard LOBPCG solver without shift-and-invert.
"""
function _solve_lobpcg_standard(solver::Solver, k, method::LOBPCGMethod; bands=1:10)
    # Build dense matrices
    LHS, RHS = build_matrices(solver, k)
    dim = size(LHS, 1)

    # Number of eigenvalues to compute
    nev = _nev(bands, dim)

    # Ensure matrices are Hermitian (required by LOBPCG)
    A = Hermitian(LHS)
    B = Hermitian(RHS)

    # Initial guess: random orthonormal vectors
    X0 = randn(ComplexF64, dim, nev)
    X0, _ = qr(X0)
    X0 = Matrix(X0)

    # Solve using LOBPCG
    # IterativeSolvers.lobpcg solves A x = λ B x for smallest eigenvalues
    results = IterativeSolvers.lobpcg(
        A, B, false, X0; tol=method.tol, maxiter=method.maxiter
    )

    # Extract eigenvalues and eigenvectors
    λ_vals = results.λ
    eigvecs = results.X

    # Convert to frequencies: ω = √λ
    ω² = real.(λ_vals)
    ω² = max.(ω², 0.0)  # Handle small negative values (numerical noise)
    frequencies = sqrt.(ω²)

    # Sort by frequency (should already be sorted, but ensure)
    perm = sortperm(frequencies)
    frequencies = frequencies[perm]
    eigenvectors = eigvecs[:, perm]

    # Select requested bands
    return _select_bands(frequencies, eigenvectors, bands)
end

"""
    _solve_lobpcg_shifted(solver, k, method; bands)

Shift-and-invert LOBPCG solver.

Uses spectral transformation to find eigenvalues above shift σ:
- Original problem: `A x = λ B x`
- Transformed: `(A - σB)⁻¹ B x = μ x` where `μ = 1/(λ - σ)`

For λ > σ: μ > 0 (positive, LOBPCG finds largest)
For λ < σ: μ < 0 (negative, filtered out)
For λ ≈ σ: |μ| → ∞ (very large, avoided by appropriate σ choice)

This is particularly useful for 3D H-field formulation where spurious
longitudinal modes exist at λ ≈ 0. Setting σ = 0.01 skips these modes.
"""
function _solve_lobpcg_shifted(solver::Solver, k, method::LOBPCGMethod; bands=1:10)
    σ = method.shift

    # Build dense matrices
    LHS, RHS = build_matrices(solver, k)
    dim = size(LHS, 1)

    # Number of eigenvalues to compute (request extra to filter)
    nev = _nev(bands, dim)
    nev_request = min(nev + 5, dim)  # Request a few extra

    # Compute shifted matrix: (A - σB)
    shifted = LHS - σ * RHS

    # LU factorization for solving (A - σB)⁻¹ B x
    shifted_fact = lu(shifted)

    # Transformed matrices for LOBPCG:
    # We want to solve: (A - σB)⁻¹ B x = μ x
    # This is equivalent to: B x = μ (A - σB) x
    # Or in standard form with A' = B, B' = (A - σB):
    # A' x = μ B' x
    # But LOBPCG needs positive definite B', which (A - σB) is NOT in general.
    #
    # Alternative approach: Convert to standard eigenvalue problem
    # (A - σB)⁻¹ B x = μ x
    # Let y = B^{1/2} x, then we need B^{1/2} (A - σB)⁻¹ B^{1/2} y = μ y
    # This requires computing B^{1/2} which can be expensive.
    #
    # Simpler approach: Use dense eigensolver on the transformed operator
    # and then convert back. This loses the LOBPCG advantage but works correctly.

    # Build the transformed operator as a dense matrix
    # T = (A - σB)⁻¹ B
    T = shifted_fact \ Matrix(RHS)

    # Make it Hermitian (should be theoretically, numerical errors may break it)
    # Note: T is NOT Hermitian in general. For LOBPCG we need Hermitian.
    # Instead, use eigen on T directly.

    # Standard eigenvalue decomposition
    μ_all, V_all = eigen(T)

    # Filter: keep only positive μ (corresponding to λ > σ)
    positive_mask = real.(μ_all) .> 0
    positive_indices = findall(positive_mask)

    if isempty(positive_indices)
        @warn "LOBPCG (shifted): No eigenvalues found above shift σ = $σ. Try a smaller shift."
        # Fall back to standard LOBPCG
        return _solve_lobpcg_standard(solver, k, method; bands=bands)
    end

    μ_pos = μ_all[positive_indices]
    V_pos = V_all[:, positive_indices]

    # Convert μ back to λ: λ = σ + 1/μ
    λ_vals = σ .+ 1.0 ./ real.(μ_pos)

    # Sort by λ (smallest first)
    perm = sortperm(λ_vals)
    λ_sorted = λ_vals[perm]
    V_sorted = V_pos[:, perm]

    # Convert to frequencies: ω = √λ
    ω² = max.(λ_sorted, 0.0)
    frequencies = sqrt.(ω²)

    # Select requested bands (uses multiple dispatch)
    nev_available = length(frequencies)
    if nev_available < nev
        @warn "LOBPCG (shifted): Only $nev_available eigenvalues found above shift σ = $σ"
    end

    freqs, vecs, warn_missing = _select_available_bands(
        frequencies, V_sorted, bands, nev_available
    )
    warn_missing && @warn "LOBPCG (shifted): Some requested bands not available"
    return freqs, vecs
end

# Dispatch for each dimension
function solve_impl(solver::Solver{Dim1}, k::Real, method::LOBPCGMethod; bands=1:10)
    _solve_lobpcg(solver, k, method; bands=bands)
end

function solve_impl(
    solver::Solver{Dim2}, k::Vector{<:Real}, method::LOBPCGMethod; bands=1:10
)
    _solve_lobpcg(solver, k, method; bands=bands)
end

function solve_impl(
    solver::Solver{Dim3}, k::Vector{<:Real}, method::LOBPCGMethod; bands=1:10
)
    _solve_lobpcg(solver, k, method; bands=bands)
end

# ============================================================================
# RSCG method implementation (placeholder for Phase 4a)
# ============================================================================

function solve_impl(solver::Solver, k, ::BasicRSCG; bands=1:10)
    error(
        "BasicRSCG is designed for Green's function / DOS / LDOS calculations, " *
        "not for eigenvalue problems. Use DenseMethod() for band structure, " *
        "or use compute_dos() / compute_ldos() with BasicRSCG solver.",
    )
end

# ============================================================================
# Fallback for unsupported methods
# ============================================================================

function solve_impl(solver::Solver, k, method::SolverMethod; bands=1:10)
    error(
        "Unsupported solver method: $(typeof(method)). " *
        "Supported methods: DenseMethod(), KrylovKitMethod(), LOBPCGMethod().",
    )
end

# ============================================================================
# Group velocity
# ============================================================================

"""
    group_velocity(solver::Solver, k; bands=1:10, δk=1e-5)

Compute group velocity v_g = ∂ω/∂k at wave vector k.

# Arguments
- `solver`: The solver
- `k`: Wave vector
- `bands`: Which bands to compute (default: 1:10)
- `δk`: Finite difference step size (default: 1e-5)

# Returns
- `v_g`: Group velocity vectors for each band
  - 2D: Vector of SVector{2} (v_gx, v_gy)
  - 1D: Vector of Float64
"""
function group_velocity(solver::Solver, k; bands=1:10, δk=1e-5)
    group_velocity_impl(solver, k, solver.method; bands=bands, δk=δk)
end

# ============================================================================
# Group velocity - Dense method (2D, finite difference)
# ============================================================================

function group_velocity_impl(
    solver::Solver{Dim2}, k::Vector{<:Real}, ::DenseMethod; bands=1:10, δk=1e-5
)
    # Central difference: ∂ω/∂k_i ≈ (ω(k+δk_i) - ω(k-δk_i)) / (2δk)

    # x-direction
    k_xp = k + [δk, 0.0]
    k_xm = k - [δk, 0.0]
    ω_xp, _ = solve(solver, k_xp; bands=bands)
    ω_xm, _ = solve(solver, k_xm; bands=bands)
    v_gx = (ω_xp - ω_xm) / (2δk)

    # y-direction
    k_yp = k + [0.0, δk]
    k_ym = k - [0.0, δk]
    ω_yp, _ = solve(solver, k_yp; bands=bands)
    ω_ym, _ = solve(solver, k_ym; bands=bands)
    v_gy = (ω_yp - ω_ym) / (2δk)

    # Return as vector of SVector{2}
    return [SVector(v_gx[i], v_gy[i]) for i in eachindex(v_gx)]
end

# Convenience methods for 2D
function group_velocity_impl(
    solver::Solver{Dim2}, k::SVector{2}, method::DenseMethod; bands=1:10, δk=1e-5
)
    group_velocity_impl(solver, Vector(k), method; bands=bands, δk=δk)
end

function group_velocity_impl(
    solver::Solver{Dim2}, k::Tuple{Real,Real}, method::DenseMethod; bands=1:10, δk=1e-5
)
    group_velocity_impl(solver, [Float64(k[1]), Float64(k[2])], method; bands=bands, δk=δk)
end

# ============================================================================
# Group velocity - Dense method (1D, finite difference)
# ============================================================================

function group_velocity_impl(
    solver::Solver{Dim1}, k::Real, ::DenseMethod; bands=1:10, δk=1e-5
)
    # Central difference: ∂ω/∂k ≈ (ω(k+δk) - ω(k-δk)) / (2δk)
    ω_p, _ = solve(solver, k + δk; bands=bands)
    ω_m, _ = solve(solver, k - δk; bands=bands)
    return (ω_p - ω_m) / (2δk)
end

# ============================================================================
# Group velocity - RSCG method (placeholder)
# ============================================================================

function group_velocity_impl(solver::Solver, k, ::BasicRSCG; bands=1:10, δk=1e-5)
    error("group_velocity with BasicRSCG is not yet implemented. Use DenseMethod().")
end

# ============================================================================
# Group velocity - Fallback
# ============================================================================

function group_velocity_impl(solver::Solver, k, method::Any; bands=1:10, δk=1e-5)
    error("Unsupported solver method for group_velocity: $(typeof(method)).")
end

# ============================================================================
# Error fallbacks for invalid argument types
# ============================================================================

"""
    Solver(wave, geometry, resolution)

Construct a solver for photonic/phononic crystal eigenvalue problems.

See concrete method signatures for detailed documentation and keyword arguments.
"""
function Solver(wave::Any, geometry::Any, resolution::Any; kwargs...)
    error(
        "Solver: expected (wave::WaveType, geometry::Geometry, resolution::Tuple{Int,...}), " *
        "got ($(typeof(wave)), $(typeof(geometry)), $(typeof(resolution)))",
    )
end

"""
    solve_at_k(solver, k, method; kwargs...)

Solve eigenvalue problem at a single k-point.

See concrete method signatures for detailed documentation and keyword arguments.
"""
function solve_at_k(solver::Any, k::Any, method::Any; kwargs...)
    error(
        "solve_at_k: expected (solver::Solver, k, method::SolverMethod), " *
        "got ($(typeof(solver)), $(typeof(k)), $(typeof(method)))",
    )
end

"""
    matrix_dimension(solver)

Return the matrix dimension for the eigenvalue problem.

See concrete method signatures for detailed documentation.
"""
function matrix_dimension(solver::Any)
    error("matrix_dimension: expected solver::Solver, got $(typeof(solver))")
end

"""
    solve_at_k_with_vectors(solver, k, method; kwargs...)

Solve eigenvalue problem and return both frequencies and eigenvectors.

See concrete method signatures for detailed documentation and keyword arguments.
"""
function solve_at_k_with_vectors(solver::Any, k::Any, method::Any; kwargs...)
    error(
        "solve_at_k_with_vectors: expected (solver::AbstractSolver, k, method::SolverMethod), " *
        "got ($(typeof(solver)), $(typeof(k)), $(typeof(method)))",
    )
end

"""
    build_matrices(solver, k)

Build the LHS and RHS matrices for the generalized eigenvalue problem.

See concrete method signatures for detailed documentation.
"""
function build_matrices(solver::Any, k::Any)
    error(
        "build_matrices: expected (solver::AbstractSolver, k), " *
        "got ($(typeof(solver)), $(typeof(k)))",
    )
end

"""
    get_weight_matrix(solver)

Return the weight matrix W for computing inner products of eigenvectors.

The weight matrix satisfies: `vecs' * W * vecs ≈ I` for orthonormal eigenvectors.
This is the RHS matrix of the generalized eigenvalue problem, which defines
the inner product in the function space.

See concrete method signatures for detailed documentation.
"""
function get_weight_matrix(solver::Any)
    error("get_weight_matrix: expected solver::AbstractSolver, got $(typeof(solver))")
end
