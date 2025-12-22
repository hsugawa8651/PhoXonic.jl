# Last-Modified: 2025-12-12T22:45:45+09:00

#=
Band structure computation and analysis.
=#

"""
    BandStructure

Results of a band structure calculation.

# Fields
- `kpoints`: K-points used (as vector of SVectors)
- `distances`: Cumulative distances along the path
- `frequencies`: Matrix of frequencies (nk × nbands)
- `labels`: High-symmetry point labels (index, label)
"""
struct BandStructure{D}
    kpoints::Vector{SVector{D,Float64}}
    distances::Vector{Float64}
    frequencies::Matrix{Float64}
    labels::Vector{Tuple{Int,String}}
end

"""
    compute_bands(solver::Solver, kpath; bands=1:10, verbose=false)

Compute band structure along a k-point path.

When the solver uses LOBPCGMethod with warm_start=true, automatically:
1. Solves first k-point with Dense (if first_dense=true)
2. Uses previous eigenvectors as initial guess for subsequent k-points
3. Applies matrix scaling (if scale=true)

This can achieve up to 38x speedup for large problems.

# Arguments
- `solver`: The solver
- `kpath`: K-point path (SimpleKPath, KPathInterpolant, or Vector)
- `bands`: Which bands to compute (default: 1:10)
- `verbose`: Print progress (default: false)

# Returns
A `BandStructure` object containing frequencies for each k-point and band.
"""
function compute_bands(solver::Solver{Dim2}, kpath::SimpleKPath{2};
                       bands=1:10, verbose::Bool=false)
    nk = length(kpath)
    nbands = length(bands)
    frequencies = zeros(Float64, nk, nbands)

    # Check if warm start is enabled
    method = solver.method
    use_warm_start = method isa LOBPCGMethod && method.warm_start

    if use_warm_start
        _compute_bands_warmstart!(frequencies, solver, kpath.points, bands, verbose)
    else
        for (ik, k) in enumerate(kpath.points)
            if verbose && ik % 10 == 0
                println("Computing k-point $ik / $nk")
            end
            ω, _ = solve(solver, collect(k); bands=bands)
            frequencies[ik, :] = ω
        end
    end

    BandStructure{2}(kpath.points, kpath.distances, frequencies, kpath.labels)
end

"""
    compute_bands(solver::Solver, kpoints::Vector; bands=1:10)

Compute band structure at specified k-points (without path structure).
"""
function compute_bands(solver::Solver{Dim2}, kpoints::Vector{<:AbstractVector};
                       bands=1:10, verbose::Bool=false)
    nk = length(kpoints)
    nbands = length(bands)
    frequencies = zeros(Float64, nk, nbands)

    for (ik, k) in enumerate(kpoints)
        if verbose && ik % 10 == 0
            println("Computing k-point $ik / $nk")
        end
        ω, _ = solve(solver, k; bands=bands)
        frequencies[ik, :] = ω
    end

    # Create simple distance array
    points = [SVector{2}(k...) for k in kpoints]
    distances = Float64[0.0]
    for i in 2:nk
        push!(distances, distances[end] + norm(points[i] - points[i-1]))
    end

    BandStructure{2}(points, distances, frequencies, Tuple{Int,String}[])
end

"""
    compute_bands(solver::Solver, kpi::KPathInterpolant; bands=1:10, verbose=false)

Compute band structure using Brillouin.jl KPathInterpolant.
"""
function compute_bands(solver::Solver{Dim2}, kpi::KPathInterpolant;
                       bands=1:10, verbose::Bool=false)
    # Extract k-points from interpolant (only kx, ky for 2D)
    kpoints_3d = collect(kpi)
    kpoints = [SVector{2}(k[1], k[2]) for k in kpoints_3d]

    nk = length(kpoints)
    nbands = length(bands)
    frequencies = zeros(Float64, nk, nbands)

    for (ik, k) in enumerate(kpoints)
        if verbose && ik % 10 == 0
            println("Computing k-point $ik / $nk")
        end
        ω, _ = solve(solver, collect(k); bands=bands)
        frequencies[ik, :] = ω
    end

    # Compute distances
    distances = Float64[0.0]
    for i in 2:nk
        push!(distances, distances[end] + norm(kpoints[i] - kpoints[i-1]))
    end

    # Extract labels from KPathInterpolant if available
    labels = Tuple{Int,String}[]

    BandStructure{2}(kpoints, distances, frequencies, labels)
end

# ============================================================================
# 3D Band structure computation
# ============================================================================

"""
    compute_bands(solver::Solver{Dim3}, kpath::SimpleKPath{3}; bands=1:10, verbose=false)

Compute 3D band structure along a k-point path.

# Arguments
- `solver`: 3D solver (FullVectorEM or FullElastic)
- `kpath`: K-point path from `simple_kpath_cubic`, `simple_kpath_fcc`, etc.
- `bands`: Which bands to compute (default: 1:10)
- `verbose`: Print progress (default: false)

# Returns
A `BandStructure{3}` object containing frequencies for each k-point and band.

# Example
```julia
lat = fcc_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Sphere([0,0,0], 0.25), Dielectric(12.0))])
solver = Solver(FullVectorEM(), geo, (12,12,12), KrylovKitMethod(shift=0.01); cutoff=3)
kpath = simple_kpath_fcc(a=1.0, npoints=20)
bands = compute_bands(solver, kpath; bands=1:6)
```

# Note
At Γ point (k=0), the lowest transverse modes also have ω→0, which can cause
anomalous values. Consider using a small k offset or skipping the Γ point.
"""
function compute_bands(solver::Solver{Dim3}, kpath::SimpleKPath{3};
                       bands=1:10, verbose::Bool=false)
    nk = length(kpath)
    nbands = length(bands)
    frequencies = zeros(Float64, nk, nbands)

    for (ik, k) in enumerate(kpath.points)
        if verbose && ik % 10 == 0
            println("Computing k-point $ik / $nk")
        end
        ω, _ = solve(solver, collect(k); bands=bands)
        frequencies[ik, :] = ω
    end

    BandStructure{3}(kpath.points, kpath.distances, frequencies, kpath.labels)
end

"""
    compute_bands(solver::Solver{Dim3}, kpoints::Vector; bands=1:10, verbose=false)

Compute 3D band structure at specified k-points (without path structure).

# Arguments
- `solver`: 3D solver (FullVectorEM or FullElastic)
- `kpoints`: Vector of k-points as 3-element vectors
- `bands`: Which bands to compute (default: 1:10)
- `verbose`: Print progress (default: false)

# Returns
A `BandStructure{3}` object. Labels will be empty since no path structure is provided.

# Example
```julia
kpoints = [[0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.5, 0.5]]
bands = compute_bands(solver, kpoints; bands=1:6)
```
"""
function compute_bands(solver::Solver{Dim3}, kpoints::Vector{<:AbstractVector};
                       bands=1:10, verbose::Bool=false)
    nk = length(kpoints)
    nbands = length(bands)
    frequencies = zeros(Float64, nk, nbands)

    for (ik, k) in enumerate(kpoints)
        if verbose && ik % 10 == 0
            println("Computing k-point $ik / $nk")
        end
        ω, _ = solve(solver, collect(k); bands=bands)
        frequencies[ik, :] = ω
    end

    # Create simple distance array
    points = [SVector{3}(k...) for k in kpoints]
    distances = Float64[0.0]
    for i in 2:nk
        push!(distances, distances[end] + norm(points[i] - points[i-1]))
    end

    BandStructure{3}(points, distances, frequencies, Tuple{Int,String}[])
end

"""
    compute_bands(solver::Solver{Dim3}, kpi::KPathInterpolant; bands=1:10, verbose=false)

Compute 3D band structure using Brillouin.jl KPathInterpolant.

# Arguments
- `solver`: 3D solver (FullVectorEM or FullElastic)
- `kpi`: K-path interpolant from `kpath_cubic`, `kpath_fcc`, `kpath_bcc`, etc.
- `bands`: Which bands to compute (default: 1:10)
- `verbose`: Print progress (default: false)

# Returns
A `BandStructure{3}` object.

# Example
```julia
lat = cubic_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Sphere([0,0,0], 0.3), Dielectric(12.0))])
solver = Solver(FullVectorEM(), geo, (12,12,12), KrylovKitMethod(shift=0.01); cutoff=3)
kpath = kpath_cubic(a=1.0, N=50)  # Brillouin.jl based
bands = compute_bands(solver, kpath; bands=1:6)
```
"""
function compute_bands(solver::Solver{Dim3}, kpi::KPathInterpolant;
                       bands=1:10, verbose::Bool=false)
    kpoints_raw = collect(kpi)
    kpoints = [SVector{3}(k[1], k[2], k[3]) for k in kpoints_raw]

    nk = length(kpoints)
    nbands = length(bands)
    frequencies = zeros(Float64, nk, nbands)

    for (ik, k) in enumerate(kpoints)
        if verbose && ik % 10 == 0
            println("Computing k-point $ik / $nk")
        end
        ω, _ = solve(solver, collect(k); bands=bands)
        frequencies[ik, :] = ω
    end

    # Compute distances
    distances = Float64[0.0]
    for i in 2:nk
        push!(distances, distances[end] + norm(kpoints[i] - kpoints[i-1]))
    end

    # Extract labels from KPathInterpolant if available
    labels = Tuple{Int,String}[]

    BandStructure{3}(kpoints, distances, frequencies, labels)
end

# ============================================================================
# Analysis functions
# ============================================================================

"""
    find_bandgap(bs::BandStructure, band1::Int, band2::Int)

Find the band gap between band1 and band2.

# Returns
Named tuple with:
- `gap`: Band gap width (0 if bands overlap)
- `gap_ratio`: Gap-to-midgap ratio
- `min_upper`: Minimum of upper band
- `max_lower`: Maximum of lower band
"""
function find_bandgap(bs::BandStructure, band1::Int, band2::Int)
    max_lower = maximum(bs.frequencies[:, band1])
    min_upper = minimum(bs.frequencies[:, band2])

    gap = max(0.0, min_upper - max_lower)
    midgap = (min_upper + max_lower) / 2
    gap_ratio = midgap > 0 ? gap / midgap : 0.0

    (gap=gap, gap_ratio=gap_ratio, min_upper=min_upper, max_lower=max_lower)
end

"""
    find_all_gaps(bs::BandStructure; threshold=0.0)

Find all band gaps in the band structure.

# Returns
Vector of named tuples with gap information.
"""
function find_all_gaps(bs::BandStructure; threshold::Float64=0.0)
    nbands = size(bs.frequencies, 2)
    gaps = []

    for i in 1:(nbands-1)
        result = find_bandgap(bs, i, i + 1)
        if result.gap > threshold
            push!(gaps, (bands=(i, i + 1), result...))
        end
    end

    return gaps
end

# ============================================================================
# Accessors
# ============================================================================

"""
    frequencies(bs::BandStructure)

Return the frequency matrix.
"""
frequencies(bs::BandStructure) = bs.frequencies

"""
    distances(bs::BandStructure)

Return the cumulative distances along the k-path.
"""
distances(bs::BandStructure) = bs.distances

"""
    labels(bs::BandStructure)

Return the high-symmetry point labels.
"""
labels(bs::BandStructure) = bs.labels

"""
    nbands(bs::BandStructure)

Return the number of bands.
"""
nbands(bs::BandStructure) = size(bs.frequencies, 2)

"""
    nkpoints(bs::BandStructure)

Return the number of k-points.
"""
nkpoints(bs::BandStructure) = size(bs.frequencies, 1)

# ============================================================================
# Warm start helper for LOBPCG
# ============================================================================

"""
    _compute_bands_warmstart!(frequencies, solver, kpoints, bands, verbose)

Internal helper for warm start band computation with LOBPCG.
Uses previous eigenvectors as initial guess for faster convergence.
"""
function _compute_bands_warmstart!(frequencies::Matrix{Float64},
                                   solver::Solver,
                                   kpoints,
                                   bands,
                                   verbose::Bool)
    nk = length(kpoints)
    method = solver.method

    # State for warm start
    prev_eigenvectors = nothing

    for (ik, k) in enumerate(kpoints)
        if verbose && ik % 10 == 0
            println("Computing k-point $ik / $nk")
        end

        k_vec = collect(k)

        if ik == 1 && method.first_dense
            # First k-point: use Dense for accurate eigenvectors
            ω, vecs = solve_at_k(solver, k_vec, DenseMethod();
                                 bands=bands, return_eigenvectors=true)
            frequencies[ik, :] = ω
            prev_eigenvectors = vecs
        else
            # Subsequent k-points: use LOBPCG with warm start
            ω, vecs = solve_at_k(solver, k_vec, method;
                                 bands=bands,
                                 X0=prev_eigenvectors,
                                 return_eigenvectors=true)
            frequencies[ik, :] = ω
            prev_eigenvectors = vecs
        end
    end
end

# ============================================================================
# Error fallback for invalid argument types
# ============================================================================

"""
    compute_bands(solver, kpath; kwargs...)

Compute band structure along a k-path.

See concrete method signatures for detailed documentation and keyword arguments.
"""
function compute_bands(solver::Any, kpath::Any; kwargs...)
    error("compute_bands: expected (solver::Solver, kpath), " *
          "got ($(typeof(solver)), $(typeof(kpath)))")
end
