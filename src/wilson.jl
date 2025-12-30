# Wilson loop and topological invariants for photonic/phononic crystals

"""
Wilson loop module for computing topological invariants.

# Exports
- `compute_zak_phase`: 1D Zak phase calculation
- `compute_wilson_spectrum`: 2D Wilson loop spectrum
- `winding_number`: Extract winding number from Wilson spectrum
"""

# Types

"""
    ZakPhaseResult

Result of 1D Zak phase calculation.

# Fields
- `phases::Vector{Float64}`: Zak phase for each band, in [-π, π]
- `bands::UnitRange{Int}`: Which bands were computed
"""
struct ZakPhaseResult
    phases::Vector{Float64}
    bands::UnitRange{Int}
end

# Core functions

"""
    overlap_matrix(a::AbstractMatrix, b::AbstractMatrix, W::AbstractMatrix)

Compute the weighted overlap matrix between two sets of eigenvectors.

The overlap matrix M[i,j] = aᵢ†Wbⱼ represents the inner product between
column i of `a` and column j of `b` with weight matrix `W`.

# Arguments
- `a`: First set of eigenvectors (n_basis × n_bands_a)
- `b`: Second set of eigenvectors (n_basis × n_bands_b)
- `W`: Weight matrix (n_basis × n_basis), typically from material parameters

# Returns
- `M`: Overlap matrix (n_bands_a × n_bands_b)
"""
function overlap_matrix(
    a::AbstractMatrix{<:Complex},
    b::AbstractMatrix{<:Complex},
    W::AbstractMatrix
)
    return a' * W * b
end

"""
    unitary_approx(M::AbstractMatrix)

Compute the unitary matrix closest to M using SVD.

Given a matrix M = UΣV†, returns UV† which is the unitary matrix
that minimizes the Frobenius norm ‖M - U‖_F.

This is used in Wilson loop calculations to extract the gauge-invariant
phase information from overlap matrices.

# Arguments
- `M`: Input matrix (n × n)

# Returns
- `U`: Unitary approximation (n × n) such that U†U = UU† = I
"""
function unitary_approx(M::AbstractMatrix{<:Complex})
    F = svd(M)
    return F.U * F.Vt
end

"""
    wilson_phases(W::AbstractMatrix)

Extract the phases from the eigenvalues of a Wilson matrix.

For a unitary Wilson matrix W, the eigenvalues lie on the unit circle
and can be written as exp(iθ). This function returns the phases θ
wrapped to the range [-π, π].

# Arguments
- `W`: Wilson matrix (n × n), typically unitary

# Returns
- `phases`: Vector of phases in [-π, π], one per eigenvalue
"""
function wilson_phases(W::AbstractMatrix{<:Complex})
    λ = eigvals(W)
    return angle.(λ)
end

"""
    wilson_matrix(spaces::Vector{<:AbstractMatrix}, W::AbstractMatrix)

Compute the Wilson matrix from a sequence of eigenvector spaces.

The Wilson matrix is the product of overlap matrices along a closed loop in k-space:
    W = U₁₂ · U₂₃ · ... · Uₙ₁
where Uᵢⱼ = unitary_approx(overlap_matrix(vᵢ, vⱼ, W)) is the unitary approximation
of the overlap between eigenvectors at adjacent k-points.

For a closed loop, the last point connects back to the first point.

# Arguments
- `spaces`: Vector of eigenvector matrices, one per k-point (n_basis × n_bands each)
- `W`: Weight matrix for the inner product (n_basis × n_basis)

# Returns
- `W_wilson`: Wilson matrix (n_bands × n_bands), unitary
"""
function wilson_matrix(
    spaces::Vector{<:AbstractMatrix{<:Complex}},
    W::AbstractMatrix
)
    n_k = length(spaces)
    n_bands = size(spaces[1], 2)

    # For a single k-point, return identity
    if n_k <= 1
        return Matrix{ComplexF64}(I, n_bands, n_bands)
    end

    # Compute product of unitary overlap matrices
    # Start with overlap between first and second k-point
    W_wilson = unitary_approx(overlap_matrix(spaces[1], spaces[2], W))

    # Multiply by subsequent overlaps
    for i in 2:(n_k-1)
        M = overlap_matrix(spaces[i], spaces[i+1], W)
        W_wilson = W_wilson * unitary_approx(M)
    end

    # Close the loop: last k-point back to first
    M_close = overlap_matrix(spaces[n_k], spaces[1], W)
    W_wilson = W_wilson * unitary_approx(M_close)

    return W_wilson
end

# High-level API

"""
    compute_zak_phase(solver::Solver{Dim1}, bands; n_k=50)

Compute the Zak phase for a 1D system.

The Zak phase is the Berry phase acquired by a Bloch state when k traverses
the full Brillouin zone from 0 to 2π/a. For systems with inversion symmetry,
the Zak phase is quantized to 0 or π.

# Arguments
- `solver`: 1D solver (Photonic1D or Longitudinal1D)
- `bands`: Which bands to compute (e.g., `1:3`)
- `n_k`: Number of k-points for integration (default: 50)

# Returns
- `ZakPhaseResult`: Contains phases for each band

# Example
```julia
lat = lattice_1d(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Segment(0.2, 0.8), Dielectric(4.0))])
solver = Solver(Photonic1D(), geo, 64; cutoff=10)
result = compute_zak_phase(solver, 1:2)
println(result.phases)  # Zak phases for bands 1 and 2
```
"""
function compute_zak_phase(
    solver::Solver{Dim1},
    bands::UnitRange{Int};
    n_k::Int = 50
)
    # Get weight matrix for inner product
    W = get_weight_matrix(solver)

    # Sample k-points from 0 to 1 (in units of reciprocal lattice vector)
    # We need n_k+1 points but the last one wraps back, so we use n_k points
    # for the Wilson loop product
    k_points = range(0.0, 1.0 - 1.0 / n_k, length = n_k)

    # Collect eigenvectors at each k-point
    n_bands = length(bands)
    spaces = Vector{Matrix{ComplexF64}}(undef, n_k)

    max_band = maximum(bands)
    for (i, k) in enumerate(k_points)
        # Solve at this k-point
        _, vectors = solve_at_k_with_vectors(solver, k, solver.method; bands = 1:max_band)

        # Extract the bands we want
        spaces[i] = vectors[:, bands]
    end

    # Compute Wilson matrix
    W_wilson = wilson_matrix(spaces, W)

    # Extract phases
    phases = wilson_phases(W_wilson)

    # Sort phases to match band ordering
    # (eigvals returns in arbitrary order, so we sort by magnitude)
    sort!(phases)

    return ZakPhaseResult(phases, bands)
end
