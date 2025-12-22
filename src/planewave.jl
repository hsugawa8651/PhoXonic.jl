# Last-Modified: 2025-12-09T22:36:15+09:00

#=
Plane wave basis for Fourier expansion.
=#

"""
    PlaneWaveBasis{D<:Dimension}

Plane wave basis for D-dimensional periodic structures.

# Fields
- `lattice`: The underlying lattice
- `cutoff`: Maximum index for plane waves
- `indices`: List of (p, q, ...) indices for each plane wave
- `G`: Reciprocal lattice vectors for each plane wave
- `num_pw`: Total number of plane waves
"""
struct PlaneWaveBasis{D<:Dimension}
    lattice::Lattice{D}
    cutoff::Int
    indices::Vector{NTuple{N,Int} where N}
    G::Vector{SVector{N,Float64} where N}
    num_pw::Int
end

# ============================================================================
# 1D Basis
# ============================================================================

"""
    PlaneWaveBasis(lattice::Lattice{Dim1}, cutoff::Int)

Create a 1D plane wave basis with indices from -cutoff to +cutoff.
"""
function PlaneWaveBasis(lattice::Lattice{Dim1}, cutoff::Int)
    b1 = lattice.reciprocal[1]

    indices = NTuple{1,Int}[]
    G_vectors = Vec1[]

    for p in (-cutoff):cutoff
        push!(indices, (p,))
        push!(G_vectors, p * b1)
    end

    PlaneWaveBasis{Dim1}(lattice, cutoff, indices, G_vectors, length(indices))
end

# ============================================================================
# 2D Basis
# ============================================================================

"""
    PlaneWaveBasis(lattice::Lattice{Dim2}, cutoff::Int)

Create a 2D plane wave basis with spherical cutoff.
Only includes indices (p, q) where p² + q² ≤ cutoff².
"""
function PlaneWaveBasis(lattice::Lattice{Dim2}, cutoff::Int)
    b1, b2 = lattice.reciprocal

    indices = NTuple{2,Int}[]
    G_vectors = Vec2[]

    for p in (-cutoff):cutoff
        for q in (-cutoff):cutoff
            if p^2 + q^2 <= cutoff^2
                push!(indices, (p, q))
                push!(G_vectors, p * b1 + q * b2)
            end
        end
    end

    PlaneWaveBasis{Dim2}(lattice, cutoff, indices, G_vectors, length(indices))
end

"""
    PlaneWaveBasis(lattice::Lattice{Dim2}, cutoff_p::Int, cutoff_q::Int)

Create a 2D plane wave basis with rectangular cutoff.
"""
function PlaneWaveBasis(lattice::Lattice{Dim2}, cutoff_p::Int, cutoff_q::Int)
    b1, b2 = lattice.reciprocal

    indices = NTuple{2,Int}[]
    G_vectors = Vec2[]

    for p in (-cutoff_p):cutoff_p
        for q in (-cutoff_q):cutoff_q
            push!(indices, (p, q))
            push!(G_vectors, p * b1 + q * b2)
        end
    end

    PlaneWaveBasis{Dim2}(
        lattice, max(cutoff_p, cutoff_q), indices, G_vectors, length(indices)
    )
end

# ============================================================================
# 3D Basis
# ============================================================================

"""
    PlaneWaveBasis(lattice::Lattice{Dim3}, cutoff::Int)

Create a 3D plane wave basis with spherical cutoff.
"""
function PlaneWaveBasis(lattice::Lattice{Dim3}, cutoff::Int)
    b1, b2, b3 = lattice.reciprocal

    indices = NTuple{3,Int}[]
    G_vectors = Vec3[]

    for p in (-cutoff):cutoff
        for q in (-cutoff):cutoff
            for r in (-cutoff):cutoff
                if p^2 + q^2 + r^2 <= cutoff^2
                    push!(indices, (p, q, r))
                    push!(G_vectors, p * b1 + q * b2 + r * b3)
                end
            end
        end
    end

    PlaneWaveBasis{Dim3}(lattice, cutoff, indices, G_vectors, length(indices))
end

# ============================================================================
# Utility functions
# ============================================================================

"""
    index_map(basis::PlaneWaveBasis)

Create a dictionary mapping indices to their position in the basis.
"""
function index_map(basis::PlaneWaveBasis)
    Dict(idx => i for (i, idx) in enumerate(basis.indices))
end

"""
    get_G(basis::PlaneWaveBasis, i::Int)

Get the i-th reciprocal lattice vector.
"""
get_G(basis::PlaneWaveBasis, i::Int) = basis.G[i]

"""
    num_planewaves(basis::PlaneWaveBasis)

Return the number of plane waves in the basis.
"""
num_planewaves(basis::PlaneWaveBasis) = basis.num_pw

Base.length(basis::PlaneWaveBasis) = basis.num_pw
