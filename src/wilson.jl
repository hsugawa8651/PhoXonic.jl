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

# Note
Per-band Zak phases have an inherent ±π gauge ambiguity due to endpoint
gauge dependence in the truncated PWE basis. The multi-band total
`sum(phases)` is gauge-invariant and reliable.
"""
struct ZakPhaseResult
    phases::Vector{Float64}
    bands::UnitRange{Int}
end

"""
    WilsonSpectrumResult

Result of 2D Wilson loop spectrum calculation.

# Fields
- `k_values::Vector{Float64}`: k-values along the scanning direction
- `phases::Matrix{Float64}`: Wilson phases (n_k_path × n_bands)
- `bands::UnitRange{Int}`: Which bands were computed
- `loop_direction::Symbol`: Direction of Wilson loop integration (:b1 or :b2)
"""
struct WilsonSpectrumResult
    k_values::Vector{Float64}
    phases::Matrix{Float64}
    bands::UnitRange{Int}
    loop_direction::Symbol
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

"""
    wilson_matrix_open(spaces, W)

Compute Wilson matrix for an open path (no closing step).

This is used when the path goes from k=0 to k=1 (inclusive), where k=1 is the
same point as k=0 due to Brillouin zone periodicity. The parallel transport
is computed through all consecutive k-points without an explicit closing step.

# Arguments
- `spaces`: Vector of eigenvector matrices, one per k-point (n_basis × n_bands each)
- `W`: Weight matrix for the inner product (n_basis × n_basis)

# Returns
- `W_wilson`: Wilson matrix (n_bands × n_bands), unitary
"""
function wilson_matrix_open(
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

    # Multiply by subsequent overlaps (all the way to the last point)
    for i in 2:(n_k-1)
        M = overlap_matrix(spaces[i], spaces[i+1], W)
        W_wilson = W_wilson * unitary_approx(M)
    end

    # No closing step - the path from k=0 to k=1 is complete
    # (k=1 = k=0 by periodicity, so this is effectively a closed loop)

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

!!! warning
    Per-band phases have ±π gauge ambiguity. Use `sum(result.phases)` for
    a gauge-invariant total.

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
    # Use n_k+1 points including endpoint k=1 for open Wilson line
    k_points = range(0.0, 1.0, length = n_k + 1)

    # Collect eigenvectors at each k-point
    n_bands = length(bands)
    spaces = Vector{Matrix{ComplexF64}}(undef, n_k + 1)

    max_band = maximum(bands)
    for (i, k) in enumerate(k_points)
        # Solve at this k-point
        _, vectors = solve_at_k_with_vectors(solver, k, solver.method; bands = 1:max_band)

        # Extract the bands we want
        spaces[i] = vectors[:, bands]
    end

    # Compute Wilson matrix
    W_wilson = wilson_matrix_open(spaces, W)

    # Extract phases
    phases = wilson_phases(W_wilson)

    # Sort phases to match band ordering
    # (eigvals returns in arbitrary order, so we sort by magnitude)
    sort!(phases)

    return ZakPhaseResult(phases, bands)
end

"""
    compute_wilson_spectrum(solver::Solver{Dim2}, bands; n_k_path=21, n_k_loop=50, loop_direction=:b2)

Compute the Wilson loop spectrum for a 2D system.

The Wilson spectrum shows the eigenvalues of the Wilson loop operator as a function
of momentum along a path in the Brillouin zone. The Wilson loop is computed along
the perpendicular direction at each k-point.

# Arguments
- `solver`: 2D solver (TEWave, TMWave, SHWave, or PSVWave)
- `bands`: Which bands to include in the Wilson loop (e.g., `1:2`)
- `n_k_path`: Number of k-points along the scanning path (default: 21)
- `n_k_loop`: Number of k-points for Wilson loop integration (default: 50)
- `loop_direction`: Direction of Wilson loop (:b1 or :b2, default: :b2)

# Returns
- `WilsonSpectrumResult`: Contains k_values, phases matrix, bands, and loop_direction

# Example
```julia
lat = square_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Circle([0.0, 0.0], 0.3), Dielectric(9.0))])
solver = Solver(TMWave(), geo, (32, 32); cutoff=5)
result = compute_wilson_spectrum(solver, 1:2; n_k_path=21, n_k_loop=50)
```
"""
function compute_wilson_spectrum(
    solver::Solver{Dim2},
    bands::UnitRange{Int};
    n_k_path::Int = 21,
    n_k_loop::Int = 50,
    loop_direction::Symbol = :b2
)
    # Get weight matrix for inner product
    W = get_weight_matrix(solver)

    n_bands = length(bands)
    max_band = maximum(bands)

    # k-values along the scanning path (0 to 1 in the scanning direction)
    k_scan = range(0.0, 1.0, length = n_k_path)

    # k-values for Wilson loop: 0 to 1 inclusive (n_k_loop+1 points)
    # This ensures proper parallel transport: k=0 → k=δ → ... → k=1
    # where k=1 is equivalent to k=0 due to BZ periodicity
    k_loop = range(0.0, 1.0, length = n_k_loop + 1)

    # Storage for Wilson phases at each scan point
    phases_matrix = Matrix{Float64}(undef, n_k_path, n_bands)

    # For each k-point along the scanning direction
    for (i_scan, k_s) in enumerate(k_scan)
        # Collect eigenvectors along the Wilson loop
        spaces = Vector{Matrix{ComplexF64}}(undef, n_k_loop + 1)

        for (i_loop, k_l) in enumerate(k_loop)
            # Construct k-vector based on loop direction
            if loop_direction == :b2
                k_vec = [k_s, k_l]  # scan in b1, loop in b2
            else  # :b1
                k_vec = [k_l, k_s]  # scan in b2, loop in b1
            end

            # Solve at this k-point
            _, vectors = solve_at_k_with_vectors(
                solver, k_vec, solver.method; bands = 1:max_band
            )

            # Extract the bands we want
            spaces[i_loop] = vectors[:, bands]
        end

        # Compute Wilson matrix for this loop (open loop, no explicit closing)
        W_wilson = wilson_matrix_open(spaces, W)

        # Extract phases and sort them
        phases = wilson_phases(W_wilson)
        sort!(phases)

        phases_matrix[i_scan, :] = phases
    end

    return WilsonSpectrumResult(collect(k_scan), phases_matrix, bands, loop_direction)
end

"""
    winding_number(result::WilsonSpectrumResult, band_index::Int)

Calculate the winding number of a Wilson phase band.

The winding number counts how many times the Wilson phase winds around
the Brillouin zone as k varies along the scanning path. A non-zero winding
number indicates non-trivial topology.

# Arguments
- `result`: Wilson spectrum result from `compute_wilson_spectrum`
- `band_index`: Which band (1 to n_bands) to compute winding for

# Returns
- `Int`: The winding number (positive for upward winding, negative for downward)

# Example
```julia
result = compute_wilson_spectrum(solver, 1:2)
w = winding_number(result, 1)  # Winding of first Wilson band
```
"""
function winding_number(result::WilsonSpectrumResult, band_index::Int)
    phases = result.phases[:, band_index]
    n_k = length(phases)

    # Compute total phase change, accounting for branch cuts
    total_wind = 0.0
    for i in 1:(n_k-1)
        Δφ = phases[i+1] - phases[i]

        # Unwrap: if jump is > π, it crossed a branch cut
        if Δφ > π
            Δφ -= 2π
        elseif Δφ < -π
            Δφ += 2π
        end

        total_wind += Δφ
    end

    # Winding number is total phase change divided by 2π, rounded to integer
    return round(Int, total_wind / (2π))
end
