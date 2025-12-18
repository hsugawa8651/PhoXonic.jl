# Last-Modified: 2025-12-09T23:13:57+09:00

#=
Lattice definitions for periodic structures.
=#

"""
    Lattice{D<:Dimension}

Represents a periodic lattice in D dimensions.

# Fields
- `vectors`: Primitive lattice vectors as a tuple of SVectors
- `reciprocal`: Reciprocal lattice vectors (computed automatically)
"""
struct Lattice{D<:Dimension}
    vectors::NTuple{N,SVector{N,Float64}} where N
    reciprocal::NTuple{N,SVector{N,Float64}} where N
end

# ============================================================================
# 1D Lattice
# ============================================================================

"""
    Lattice(a::Real)

Create a 1D lattice with period `a`.
"""
function Lattice(a::Real)
    v = (Vec1(a),)
    r = (Vec1(2π / a),)
    Lattice{Dim1}(v, r)
end

# ============================================================================
# 2D Lattice
# ============================================================================

"""
    Lattice(a1::Vec2, a2::Vec2)

Create a 2D lattice with primitive vectors `a1` and `a2`.
"""
function Lattice(a1::Vec2, a2::Vec2)
    # Compute reciprocal lattice vectors
    # b1 · a1 = 2π, b1 · a2 = 0
    # b2 · a1 = 0,  b2 · a2 = 2π
    V = a1[1] * a2[2] - a1[2] * a2[1]  # Area of unit cell
    b1 = (2π / V) * Vec2(a2[2], -a2[1])
    b2 = (2π / V) * Vec2(-a1[2], a1[1])
    Lattice{Dim2}((a1, a2), (b1, b2))
end

"""
    Lattice(a1::AbstractVector, a2::AbstractVector)

Create a 2D lattice from any vector types (converted to Vec2).
"""
function Lattice(a1::AbstractVector, a2::AbstractVector)
    Lattice(Vec2(a1...), Vec2(a2...))
end

# ============================================================================
# 3D Lattice
# ============================================================================

"""
    Lattice(a1::Vec3, a2::Vec3, a3::Vec3)

Create a 3D lattice with primitive vectors `a1`, `a2`, and `a3`.
"""
function Lattice(a1::Vec3, a2::Vec3, a3::Vec3)
    # Compute reciprocal lattice vectors using cross products
    # b_i = 2π (a_j × a_k) / (a_i · (a_j × a_k))
    V = dot(a1, cross(a2, a3))  # Volume of unit cell
    b1 = (2π / V) * cross(a2, a3)
    b2 = (2π / V) * cross(a3, a1)
    b3 = (2π / V) * cross(a1, a2)
    Lattice{Dim3}((a1, a2, a3), (b1, b2, b3))
end

# ============================================================================
# Convenience constructors
# ============================================================================

"""
    lattice_1d(a::Real=1.0)

Create a 1D periodic lattice with period `a`.
"""
lattice_1d(a::Real=1.0) = Lattice(Float64(a))

"""
    square_lattice(a::Real=1.0)

Create a 2D square lattice with lattice constant `a`.
"""
function square_lattice(a::Real=1.0)
    Lattice(Vec2(a, 0), Vec2(0, a))
end

"""
    hexagonal_lattice(a::Real=1.0)

Create a 2D hexagonal (triangular) lattice with lattice constant `a`.
"""
function hexagonal_lattice(a::Real=1.0)
    Lattice(Vec2(a, 0), Vec2(a / 2, a * sqrt(3) / 2))
end

"""
    cubic_lattice(a::Real=1.0)

Create a 3D simple cubic lattice with lattice constant `a`.
"""
function cubic_lattice(a::Real=1.0)
    Lattice(Vec3(a, 0, 0), Vec3(0, a, 0), Vec3(0, 0, a))
end

"""
    fcc_lattice(a::Real=1.0)

Create a 3D face-centered cubic lattice with conventional lattice constant `a`.
"""
function fcc_lattice(a::Real=1.0)
    Lattice(
        Vec3(0, a/2, a/2),
        Vec3(a/2, 0, a/2),
        Vec3(a/2, a/2, 0)
    )
end

# ============================================================================
# Accessor functions
# ============================================================================

"""
    reciprocal_vectors(lattice::Lattice)

Return the reciprocal lattice vectors.
"""
reciprocal_vectors(lattice::Lattice) = lattice.reciprocal

"""
    lattice_vectors(lattice::Lattice)

Return the primitive lattice vectors.
"""
lattice_vectors(lattice::Lattice) = lattice.vectors

"""
    cell_area(lattice::Lattice{Dim2})

Return the area of the 2D unit cell.
"""
function cell_area(lattice::Lattice{Dim2})
    a1, a2 = lattice.vectors
    abs(a1[1] * a2[2] - a1[2] * a2[1])
end

"""
    cell_volume(lattice::Lattice{Dim3})

Return the volume of the 3D unit cell.
"""
function cell_volume(lattice::Lattice{Dim3})
    a1, a2, a3 = lattice.vectors
    abs(dot(a1, cross(a2, a3)))
end
