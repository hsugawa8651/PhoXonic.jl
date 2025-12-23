# Last-Modified: 2025-12-12T22:47:30+09:00

#=
k-point path generation using Brillouin.jl.

This module provides integration with Brillouin.jl for standardized
k-path generation through the Brillouin zone.
=#

using Brillouin: KPath as BrillouinKPath, KPathInterpolant
using Brillouin: irrfbz_path, interpolate, cartesianize

# Re-export Brillouin types for convenience
export BrillouinKPath, KPathInterpolant

"""
    kpath_from_brillouin(sgnum::Int, Rs; N=100)

Create a k-path using Brillouin.jl's irrfbz_path.

# Arguments
- `sgnum`: Space group number (e.g., 1 for P1, 227 for diamond)
- `Rs`: Direct lattice vectors as a matrix or vector of vectors
- `N`: Number of interpolation points (default: 100)

# Returns
A KPathInterpolant that can be iterated over.

# Example
```julia
# Square lattice (space group 1, primitive cell)
Rs = [[1.0, 0.0], [0.0, 1.0]]
kpi = kpath_from_brillouin(1, Rs; N=100)

# Iterate over k-points
for k in kpi
    # k is an SVector
end
```
"""
function kpath_from_brillouin(sgnum::Int, Rs; N::Int=100)
    kp = irrfbz_path(sgnum, Rs)
    kpi = interpolate(kp, N)
    return cartesianize(kpi)
end

"""
    kpath_square(; a=1.0, N=100)

Create a standard k-path for a 2D square lattice: Γ → X → M → Γ

Uses Brillouin.jl with space group 1 (P1).
"""
function kpath_square(; a::Real=1.0, N::Int=100)
    Rs = [[a, 0.0], [0.0, a]]
    # For 2D, we use space group 1 with 2D lattice vectors
    # Brillouin.jl requires 3D, so we add a z-component
    Rs_3d = [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, 1.0]]
    kp = irrfbz_path(1, Rs_3d)
    kpi = interpolate(kp, N)
    return cartesianize(kpi)
end

"""
    kpath_hexagonal(; a=1.0, N=100)

Create a standard k-path for a 2D hexagonal lattice: Γ → M → K → Γ

Uses Brillouin.jl with appropriate space group.
"""
function kpath_hexagonal(; a::Real=1.0, N::Int=100)
    # Hexagonal lattice vectors
    Rs_3d = [[a, 0.0, 0.0], [a / 2, a * sqrt(3) / 2, 0.0], [0.0, 0.0, 1.0]]
    kp = irrfbz_path(1, Rs_3d)
    kpi = interpolate(kp, N)
    return cartesianize(kpi)
end

# ============================================================================
# Simple fallback k-path for cases where Brillouin.jl is not needed
# ============================================================================

"""
    SimpleKPath{D}

A simple k-path structure for basic band structure calculations.
Use this when Brillouin.jl's full functionality is not needed.
"""
struct SimpleKPath{D}
    points::Vector{SVector{D,Float64}}
    labels::Vector{Tuple{Int,String}}
    distances::Vector{Float64}
end

"""
    SimpleKPath(waypoints, labels; npoints=50)

Create a simple k-path through specified waypoints.
"""
function SimpleKPath(
    waypoints::Vector{<:AbstractVector}, labels::Vector{String}; npoints::Int=50
)
    # Input validation
    if isempty(waypoints)
        error("SimpleKPath: waypoints must not be empty")
    end
    if length(waypoints) != length(labels)
        error("SimpleKPath: waypoints and labels must have the same length")
    end
    if npoints < 1
        error("SimpleKPath: npoints must be at least 1")
    end

    D = length(waypoints[1])

    all_points = SVector{D,Float64}[]
    all_labels = Tuple{Int,String}[]
    all_distances = Float64[]

    current_dist = 0.0

    for i in 1:length(waypoints)
        if i == 1
            push!(all_points, SVector{D}(waypoints[1]...))
            push!(all_labels, (1, labels[1]))
            push!(all_distances, 0.0)
        else
            k_start = SVector{D}(waypoints[i - 1]...)
            k_end = SVector{D}(waypoints[i]...)
            segment_length = norm(k_end - k_start)

            for j in 1:npoints
                t = j / npoints
                k = k_start + t * (k_end - k_start)
                push!(all_points, k)
                push!(all_distances, current_dist + t * segment_length)
            end

            current_dist += segment_length
            push!(all_labels, (length(all_points), labels[i]))
        end
    end

    SimpleKPath{D}(all_points, all_labels, all_distances)
end

"""
    simple_kpath_square(; a=1.0, npoints=50)

Create a simple k-path for square lattice without Brillouin.jl.
Γ → X → M → Γ
"""
function simple_kpath_square(; a::Real=1.0, npoints::Int=50)
    scale = 2π / a
    Γ = [0.0, 0.0]
    X = [0.5 * scale, 0.0]
    M = [0.5 * scale, 0.5 * scale]

    SimpleKPath([Γ, X, M, Γ], ["Γ", "X", "M", "Γ"]; npoints=npoints)
end

"""
    simple_kpath_hexagonal(; a=1.0, npoints=50)

Create a simple k-path for hexagonal lattice without Brillouin.jl.
Γ → M → K → Γ
"""
function simple_kpath_hexagonal(; a::Real=1.0, npoints::Int=50)
    # Reciprocal lattice vectors for hexagonal
    b1 = (2π / a) * SVector(1.0, -1.0 / sqrt(3))
    b2 = (2π / a) * SVector(0.0, 2.0 / sqrt(3))

    Γ = [0.0, 0.0]
    M = collect(0.5 * b1)
    K = collect((1.0 / 3.0) * b1 + (1.0 / 3.0) * b2)

    SimpleKPath([Γ, M, K, Γ], ["Γ", "M", "K", "Γ"]; npoints=npoints)
end

# Iteration interface for SimpleKPath
Base.length(kp::SimpleKPath) = length(kp.points)
Base.getindex(kp::SimpleKPath, i::Int) = kp.points[i]
Base.iterate(kp::SimpleKPath) = iterate(kp.points)
Base.iterate(kp::SimpleKPath, state) = iterate(kp.points, state)

# ============================================================================
# 3D k-paths using Brillouin.jl
# ============================================================================

"""
    kpath_cubic(; a=1.0, N=100)

Create a k-path for simple cubic (SC) lattice using Brillouin.jl.

Path: Γ → X → M → Γ → R → X | R → M

Uses space group 221 (Pm-3m).

# Arguments
- `a`: Lattice constant (default: 1.0)
- `N`: Number of interpolation points (default: 100)

# Returns
A `KPathInterpolant` from Brillouin.jl that can be passed to `compute_bands`.

# Example
```julia
lat = cubic_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Sphere([0,0,0], 0.3), Dielectric(12.0))])
solver = Solver(FullVectorEM(), geo, (12,12,12), KrylovKitMethod(shift=0.01); cutoff=3)
kpath = kpath_cubic(a=1.0, N=100)
bands = compute_bands(solver, kpath; bands=1:6)
```

Reference: Setyawan & Curtarolo, Comp. Mat. Sci. 49, 299 (2010). DOI:10.1016/j.commatsci.2010.05.010
"""
function kpath_cubic(; a::Real=1.0, N::Int=100)
    Rs = [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]]
    kp = irrfbz_path(221, Rs)  # Pm-3m
    kpi = interpolate(kp, N)
    return cartesianize(kpi)
end

"""
    kpath_fcc(; a=1.0, N=100)

Create a k-path for face-centered cubic (FCC) lattice using Brillouin.jl.

Path: Γ → X → U | K → Γ → L → W → X

Uses space group 225 (Fm-3m).

# Arguments
- `a`: Conventional lattice constant (default: 1.0)
- `N`: Number of interpolation points (default: 100)

# Returns
A `KPathInterpolant` from Brillouin.jl that can be passed to `compute_bands`.

# Example
```julia
lat = fcc_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Sphere([0,0,0], 0.25), Dielectric(12.0))])
solver = Solver(FullVectorEM(), geo, (12,12,12), KrylovKitMethod(shift=0.01); cutoff=3)
kpath = kpath_fcc(a=1.0, N=100)
bands = compute_bands(solver, kpath; bands=1:6)
```

Reference: Setyawan & Curtarolo, Comp. Mat. Sci. 49, 299 (2010). DOI:10.1016/j.commatsci.2010.05.010
"""
function kpath_fcc(; a::Real=1.0, N::Int=100)
    # FCC primitive vectors
    Rs = [[0.0, a/2, a/2], [a/2, 0.0, a/2], [a/2, a/2, 0.0]]
    kp = irrfbz_path(225, Rs)  # Fm-3m
    kpi = interpolate(kp, N)
    return cartesianize(kpi)
end

"""
    kpath_bcc(; a=1.0, N=100)

Create a k-path for body-centered cubic (BCC) lattice using Brillouin.jl.

Path: Γ → H → N → Γ → P → H | P → N

Uses space group 229 (Im-3m).

# Arguments
- `a`: Conventional lattice constant (default: 1.0)
- `N`: Number of interpolation points (default: 100)

# Returns
A `KPathInterpolant` from Brillouin.jl that can be passed to `compute_bands`.

# Example
```julia
kpath = kpath_bcc(a=1.0, N=100)
bands = compute_bands(solver, kpath; bands=1:6)
```

Reference: Setyawan & Curtarolo, Comp. Mat. Sci. 49, 299 (2010). DOI:10.1016/j.commatsci.2010.05.010
"""
function kpath_bcc(; a::Real=1.0, N::Int=100)
    # BCC primitive vectors
    Rs = [[-a/2, a/2, a/2], [a/2, -a/2, a/2], [a/2, a/2, -a/2]]
    kp = irrfbz_path(229, Rs)  # Im-3m
    kpi = interpolate(kp, N)
    return cartesianize(kpi)
end

# ============================================================================
# Simple 3D k-paths (without Brillouin.jl dependency, single connected path)
# ============================================================================

"""
    simple_kpath_cubic(; a=1.0, npoints=50)

Create a simple k-path for simple cubic (SC) lattice.

Path: Γ → X → M → Γ → R → X

High-symmetry points (in units of 2π/a):
- Γ = (0, 0, 0)
- X = (1/2, 0, 0)
- M = (1/2, 1/2, 0)
- R = (1/2, 1/2, 1/2)

# Arguments
- `a`: Lattice constant (default: 1.0)
- `npoints`: Number of points per segment (default: 50)

# Returns
A `SimpleKPath{3}` object that can be passed to `compute_bands`.

# Example
```julia
lat = cubic_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Sphere([0,0,0], 0.3), Dielectric(12.0))])
solver = Solver(FullVectorEM(), geo, (12,12,12), KrylovKitMethod(shift=0.01); cutoff=3)
kpath = simple_kpath_cubic(a=1.0, npoints=30)
bands = compute_bands(solver, kpath; bands=1:6)
```
"""
function simple_kpath_cubic(; a::Real=1.0, npoints::Int=50)
    scale = 2π / a
    Γ = [0.0, 0.0, 0.0]
    X = [0.5 * scale, 0.0, 0.0]
    M = [0.5 * scale, 0.5 * scale, 0.0]
    R = [0.5 * scale, 0.5 * scale, 0.5 * scale]

    SimpleKPath([Γ, X, M, Γ, R, X], ["Γ", "X", "M", "Γ", "R", "X"]; npoints=npoints)
end

"""
    simple_kpath_fcc(; a=1.0, npoints=50)

Create a simple k-path for face-centered cubic (FCC) lattice.

Path: Γ → X → W → L → Γ → K

High-symmetry points (values shown are multiplied by 2π/a to give Cartesian coordinates):
- Γ = (0, 0, 0)
- X = (0, 1, 0) × 2π/a
- W = (1/2, 1, 0) × 2π/a
- L = (1/2, 1/2, 1/2) × 2π/a
- K = (3/4, 3/4, 0) × 2π/a

# Arguments
- `a`: Conventional lattice constant (default: 1.0)
- `npoints`: Number of points per segment (default: 50)

# Returns
A `SimpleKPath{3}` object that can be passed to `compute_bands`.

# Example
```julia
lat = fcc_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Sphere([0,0,0], 0.25), Dielectric(12.0))])
solver = Solver(FullVectorEM(), geo, (12,12,12), KrylovKitMethod(shift=0.01); cutoff=3)
kpath = simple_kpath_fcc(a=1.0, npoints=30)
bands = compute_bands(solver, kpath; bands=1:6)
```
"""
function simple_kpath_fcc(; a::Real=1.0, npoints::Int=50)
    scale = 2π / a
    Γ = [0.0, 0.0, 0.0]
    X = [1.0 * scale, 0.0, 0.0]           # (1, 0, 0) × 2π/a
    W = [1.0 * scale, 0.5 * scale, 0.0]   # (1, 1/2, 0) × 2π/a
    L = [0.5 * scale, 0.5 * scale, 0.5 * scale]  # (1/2, 1/2, 1/2) × 2π/a
    K = [0.75 * scale, 0.75 * scale, 0.0]        # (3/4, 3/4, 0) × 2π/a

    SimpleKPath([Γ, X, W, L, Γ, K], ["Γ", "X", "W", "L", "Γ", "K"]; npoints=npoints)
end

"""
    simple_kpath_bcc(; a=1.0, npoints=50)

Create a simple k-path for body-centered cubic (BCC) lattice.

Path: Γ → H → N → Γ → P → H

High-symmetry points (values shown are multiplied by 2π/a to give Cartesian coordinates):
- Γ = (0, 0, 0)
- H = (0, 0, 1) × 2π/a
- N = (0, 1/2, 1/2) × 2π/a
- P = (1/4, 1/4, 1/4) × 2π/a

# Arguments
- `a`: Conventional lattice constant (default: 1.0)
- `npoints`: Number of points per segment (default: 50)

# Returns
A `SimpleKPath{3}` object that can be passed to `compute_bands`.

# Example
```julia
kpath = simple_kpath_bcc(a=1.0, npoints=30)
bands = compute_bands(solver, kpath; bands=1:6)
```
"""
function simple_kpath_bcc(; a::Real=1.0, npoints::Int=50)
    scale = 2π / a
    Γ = [0.0, 0.0, 0.0]
    H = [0.0, 0.0, 1.0 * scale]
    N = [0.0, 0.5 * scale, 0.5 * scale]
    P = [0.25 * scale, 0.25 * scale, 0.25 * scale]

    SimpleKPath([Γ, H, N, Γ, P, H], ["Γ", "H", "N", "Γ", "P", "H"]; npoints=npoints)
end
