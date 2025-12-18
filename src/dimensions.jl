# Last-Modified: 2025-12-09T22:32:21+09:00

#=
Dimension type system for 1D, 2D, and 3D structures.
=#

"""
Abstract type representing spatial dimensions.
"""
abstract type Dimension end

"""
    Dim1

One-dimensional structures (e.g., Bragg mirrors, superlattices).
"""
struct Dim1 <: Dimension end

"""
    Dim2

Two-dimensional structures (e.g., photonic crystal slabs, 2D phononic crystals).
"""
struct Dim2 <: Dimension end

"""
    Dim3

Three-dimensional structures (e.g., 3D photonic crystals, bulk phononic crystals).
"""
struct Dim3 <: Dimension end

# Type aliases for static vectors
const Vec1 = SVector{1,Float64}
const Vec2 = SVector{2,Float64}
const Vec3 = SVector{3,Float64}

"""
    VecType(::Type{D}) where D <: Dimension

Return the appropriate vector type for dimension D.
"""
VecType(::Type{Dim1}) = Vec1
VecType(::Type{Dim2}) = Vec2
VecType(::Type{Dim3}) = Vec3

"""
    ndims(::Type{D}) where D <: Dimension

Return the number of spatial dimensions.
"""
Base.ndims(::Type{Dim1}) = 1
Base.ndims(::Type{Dim2}) = 2
Base.ndims(::Type{Dim3}) = 3
