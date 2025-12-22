# Last-Modified: 2025-12-15T00:38:00+09:00

#=
Geometry definition for periodic structures.
=#

"""
    Geometry{D<:Dimension, M<:Material}

Represents the geometry of a periodic crystal structure.

# Fields
- `lattice`: The periodic lattice
- `background`: Background material
- `inclusions`: List of (shape, material) pairs

# Examples
```julia
lattice = square_lattice(1.0)
air = Dielectric(1.0)
rod = Dielectric(8.9)
geo = Geometry(lattice, air, [(Circle(Vec2(0,0), 0.2), rod)])
```
"""
struct Geometry{D<:Dimension,M<:Material}
    lattice::Lattice{D}
    background::M
    inclusions::Vector{Tuple{Shape{D},M}}
end

"""
    Geometry(lattice, background)

Create a homogeneous geometry with no inclusions.
"""
function Geometry(lattice::Lattice{D}, background::M) where {D<:Dimension,M<:Material}
    Geometry{D,M}(lattice, background, Tuple{Shape{D},M}[])
end

"""
    Geometry(lattice, background, inclusions)

Create a geometry with inclusions.
"""
function Geometry(lattice::Lattice{D}, background::M,
                  inclusions::Vector{<:Tuple{<:Shape{D},M2}}) where {D,M<:Material,M2<:Material}
    # Ensure all materials have compatible types
    MT = promote_type(M, M2)
    incl = Vector{Tuple{Shape{D},MT}}([(s, convert(MT, m)) for (s, m) in inclusions])
    Geometry{D,MT}(lattice, convert(MT, background), incl)
end

# Helper to extract dimension type from Shape
_shape_dim(::Shape{D}) where D = D
_dim_name(::Type{Dim1}) = "1D"
_dim_name(::Type{Dim2}) = "2D"
_dim_name(::Type{Dim3}) = "3D"

"""
    Geometry(lattice, background, inclusions) - fallback for dimension mismatch

Provides a helpful error message when shape dimensions don't match lattice dimension.
"""
function Geometry(lattice::Lattice{D}, background::M,
                  inclusions::Vector{<:Tuple{<:Shape,M2}}) where {D<:Dimension,M<:Material,M2<:Material}
    lat_dim = _dim_name(D)
    for (i, (shape, _)) in enumerate(inclusions)
        shape_dim_type = _shape_dim(shape)
        if shape_dim_type !== D
            shape_dim = _dim_name(shape_dim_type)
            shape_name = nameof(typeof(shape))
            throw(ArgumentError(
                "Dimension mismatch: cannot add $shape_dim shape ($shape_name) " *
                "to $lat_dim geometry."
            ))
        end
    end
    # If dimensions all match, use promote_type for mixed materials
    MT = promote_type(M, M2)
    typed_incl = Vector{Tuple{Shape{D},MT}}([(s, convert(MT, m)) for (s, m) in inclusions])
    Geometry{D,MT}(lattice, convert(MT, background), typed_incl)
end

"""
    get_material(geo::Geometry, point)

Return the material at a given point in real space.
Checks inclusions in order; first match wins.

Uses periodic wrapping to ensure points are correctly assigned to shapes
even when shapes are centered at lattice corners (e.g., center = (0,0)).
"""
function get_material(geo::Geometry{D}, point) where D
    for (shape, material) in geo.inclusions
        if point_in_shape_periodic(point, shape, geo.lattice)
            return material
        end
    end
    return geo.background
end

"""
    point_in_shape_periodic(point, shape, lattice)

Check if a point is inside a shape, considering periodic boundary conditions.
This wraps the point relative to the shape's center before checking.
"""
function point_in_shape_periodic end

# Default: just use the standard check (for shapes without known centers)
point_in_shape_periodic(point, shape, lattice) = point in shape

# Circle: wrap point relative to circle center
function point_in_shape_periodic(point::Vec2, c::Circle, lattice::Lattice{Dim2})
    # Get lattice vectors (for wrapping)
    a1, a2 = lattice.vectors

    # Compute displacement from circle center
    d = point - c.center

    # For non-orthogonal lattices, we need to check all nearby periodic images
    # The minimum image convention: find the closest periodic image
    r² = c.radius^2

    # Check the original position and all 8 neighboring periodic images
    for i in -1:1
        for j in -1:1
            # Displaced position: d + i*a1 + j*a2
            d_shifted = d + i * a1 + j * a2
            if dot(d_shifted, d_shifted) <= r²
                return true
            end
        end
    end

    return false
end

function point_in_shape_periodic(point::AbstractVector, c::Circle, lattice::Lattice{Dim2})
    point_in_shape_periodic(Vec2(point...), c, lattice)
end

# Rectangle: wrap point relative to rectangle center
function point_in_shape_periodic(point::Vec2, rect::Rectangle, lattice::Lattice{Dim2})
    a1, a2 = lattice.vectors

    d = point - rect.center
    hw = rect.size[1] / 2  # half-width
    hh = rect.size[2] / 2  # half-height

    # Check original and all 8 neighboring periodic images
    for i in -1:1
        for j in -1:1
            d_shifted = d + i * a1 + j * a2
            if abs(d_shifted[1]) <= hw && abs(d_shifted[2]) <= hh
                return true
            end
        end
    end

    return false
end

function point_in_shape_periodic(point::AbstractVector, r::Rectangle, lattice::Lattice{Dim2})
    point_in_shape_periodic(Vec2(point...), r, lattice)
end

# Polygon: check if point is in polygon or any periodic image
function point_in_shape_periodic(point::Vec2, poly::Polygon, lattice::Lattice{Dim2})
    a1, a2 = lattice.vectors

    # Check original and all 8 neighboring periodic images
    for i in -1:1
        for j in -1:1
            shifted_point = point - i * a1 - j * a2
            if shifted_point in poly
                return true
            end
        end
    end

    return false
end

function point_in_shape_periodic(point::AbstractVector, poly::Polygon, lattice::Lattice{Dim2})
    point_in_shape_periodic(Vec2(point...), poly, lattice)
end

# 3D Sphere: wrap point relative to sphere center
function point_in_shape_periodic(point::Vec3, s::Sphere, lattice::Lattice{Dim3})
    a1, a2, a3 = lattice.vectors
    r² = s.radius^2

    d = point - s.center

    # Check original and all 26 neighboring periodic images (3x3x3 - 1 + 1 = 27)
    for i in -1:1
        for j in -1:1
            for k in -1:1
                d_shifted = d + i * a1 + j * a2 + k * a3
                if dot(d_shifted, d_shifted) <= r²
                    return true
                end
            end
        end
    end

    return false
end

function point_in_shape_periodic(point::AbstractVector, s::Sphere, lattice::Lattice{Dim3})
    point_in_shape_periodic(Vec3(point...), s, lattice)
end

# 3D Cylinder: check all periodic images
function point_in_shape_periodic(point::Vec3, cyl::Cylinder, lattice::Lattice{Dim3})
    a1, a2, a3 = lattice.vectors

    # Check original and all 26 neighboring periodic images
    for i in -1:1
        for j in -1:1
            for k in -1:1
                shifted_point = point - i * a1 - j * a2 - k * a3
                if shifted_point in cyl
                    return true
                end
            end
        end
    end

    return false
end

function point_in_shape_periodic(point::AbstractVector, cyl::Cylinder, lattice::Lattice{Dim3})
    point_in_shape_periodic(Vec3(point...), cyl, lattice)
end

# 3D Slab: only check z-direction periodic images (infinite in xy)
function point_in_shape_periodic(point::Vec3, slab::Slab, lattice::Lattice{Dim3})
    a3 = lattice.vectors[3]

    # Check original and neighboring periodic images in z-direction
    for k in -1:1
        shifted_point = point - k * a3
        if shifted_point in slab
            return true
        end
    end

    return false
end

function point_in_shape_periodic(point::AbstractVector, slab::Slab, lattice::Lattice{Dim3})
    point_in_shape_periodic(Vec3(point...), slab, lattice)
end

# ============================================================================
# Discretization
# ============================================================================

"""
    DiscretizationMethod

Abstract type for material discretization methods.
"""
abstract type DiscretizationMethod end

"""
    SimpleGrid

Simple point-sampling discretization.
Each grid point gets the material value at that exact location.
"""
struct SimpleGrid <: DiscretizationMethod end

"""
    SubpixelAverage

Subpixel averaging for smoother material boundaries.
Uses volume fractions to average material properties.
"""
struct SubpixelAverage <: DiscretizationMethod
    samples_per_dim::Int  # Number of samples per dimension in each cell
end

SubpixelAverage() = SubpixelAverage(4)

"""
    discretize(geo::Geometry{Dim2}, resolution, property, [method])

Discretize a 2D geometry onto a grid.

# Arguments
- `geo`: The geometry to discretize
- `resolution`: Tuple (Nx, Ny) of grid points
- `property`: Property to extract (:ε, :μ, :ρ, :C11, :C12, :C44)
- `method`: Discretization method (default: SimpleGrid())

# Returns
A matrix of property values on the grid.
"""
function discretize(geo::Geometry{Dim2}, resolution::Tuple{Int,Int},
                    property::Symbol, method::DiscretizationMethod=SimpleGrid())
    discretize_impl(geo, resolution, property, method)
end

function discretize_impl(geo::Geometry{Dim2}, resolution::Tuple{Int,Int},
                         property::Symbol, ::SimpleGrid)
    Nx, Ny = resolution
    a1, a2 = geo.lattice.vectors

    result = zeros(Float64, Nx, Ny)

    for j in 1:Ny, i in 1:Nx
        # Fractional coordinates [0, 1)
        fx = (i - 1) / Nx
        fy = (j - 1) / Ny
        # Real space coordinates
        point = fx * a1 + fy * a2
        # Get material at this point
        mat = get_material(geo, point)
        # Extract requested property
        result[i, j] = get_property(mat, property)
    end

    return result
end

function discretize_impl(geo::Geometry{Dim2}, resolution::Tuple{Int,Int},
                         property::Symbol, method::SubpixelAverage)
    Nx, Ny = resolution
    a1, a2 = geo.lattice.vectors
    ns = method.samples_per_dim

    result = zeros(Float64, Nx, Ny)

    for j in 1:Ny, i in 1:Nx
        # Average over subpixel samples
        val = 0.0
        for sj in 1:ns, si in 1:ns
            fx = (i - 1 + (si - 0.5) / ns) / Nx
            fy = (j - 1 + (sj - 0.5) / ns) / Ny
            point = fx * a1 + fy * a2
            mat = get_material(geo, point)
            val += get_property(mat, property)
        end
        result[i, j] = val / (ns * ns)
    end

    return result
end

"""
    get_property(mat::Material, property::Symbol)

Extract a specific property from a material.
"""
# Error fallback for unsupported material types
get_property(mat, property::Symbol) =
    throw(ArgumentError("Unsupported material type: $(typeof(mat))"))

function get_property(mat::Dielectric, property::Symbol)
    if property == :ε
        return mat.ε
    elseif property == :μ
        return mat.μ
    elseif property == :ε_inv
        return 1.0 / mat.ε
    elseif property == :μ_inv
        return 1.0 / mat.μ
    else
        error("Unknown property $property for Dielectric")
    end
end

function get_property(mat::IsotropicElastic, property::Symbol)
    if property == :ρ
        return mat.ρ
    elseif property == :C11
        return mat.C11
    elseif property == :C12
        return mat.C12
    elseif property == :C44
        return mat.C44
    else
        error("Unknown property $property for IsotropicElastic")
    end
end

function get_property(mat::ElasticVoid, property::Symbol)
    if property == :ρ
        return mat.ρ
    elseif property == :C11
        return mat.C11
    elseif property == :C12
        return mat.C12
    elseif property == :C44
        return mat.C44
    else
        error("Unknown property $property for ElasticVoid")
    end
end

# ============================================================================
# 1D Discretization
# ============================================================================

"""
    discretize(geo::Geometry{Dim1}, resolution, property, [method])

Discretize a 1D geometry onto a grid.
"""
function discretize(geo::Geometry{Dim1}, resolution::Tuple{Int},
                    property::Symbol, method::DiscretizationMethod=SimpleGrid())
    discretize_impl(geo, resolution, property, method)
end

# Also accept Int directly for 1D
function discretize(geo::Geometry{Dim1}, resolution::Int,
                    property::Symbol, method::DiscretizationMethod=SimpleGrid())
    discretize_impl(geo, (resolution,), property, method)
end

function discretize_impl(geo::Geometry{Dim1}, resolution::Tuple{Int},
                         property::Symbol, ::SimpleGrid)
    Nx = resolution[1]
    a = geo.lattice.vectors[1][1]  # Lattice constant

    result = zeros(Float64, Nx)

    for i in 1:Nx
        # Fractional coordinate [0, 1)
        fx = (i - 1) / Nx
        # Real space coordinate
        x = fx * a
        # Get material at this point
        mat = get_material(geo, x)
        # Extract requested property
        result[i] = get_property(mat, property)
    end

    return result
end

function discretize_impl(geo::Geometry{Dim1}, resolution::Tuple{Int},
                         property::Symbol, method::SubpixelAverage)
    Nx = resolution[1]
    a = geo.lattice.vectors[1][1]
    ns = method.samples_per_dim

    result = zeros(Float64, Nx)

    for i in 1:Nx
        val = 0.0
        for si in 1:ns
            fx = (i - 1 + (si - 0.5) / ns) / Nx
            x = fx * a
            mat = get_material(geo, x)
            val += get_property(mat, property)
        end
        result[i] = val / ns
    end

    return result
end

# 1D get_material with scalar x (with periodic boundary handling)
function get_material(geo::Geometry{Dim1}, x::Real)
    for (shape, material) in geo.inclusions
        if point_in_shape_periodic_1d(x, shape, geo.lattice)
            return material
        end
    end
    return geo.background
end

"""
    point_in_shape_periodic_1d(x, shape, lattice)

Check if a 1D point is inside a shape, considering periodic boundary conditions.
"""
function point_in_shape_periodic_1d end

# Default: just use the standard check
point_in_shape_periodic_1d(x, shape, lattice) = x in shape

# Segment: check x and periodic images x ± a
function point_in_shape_periodic_1d(x::Real, seg::Segment, lattice::Lattice{Dim1})
    a = lattice.vectors[1][1]  # Lattice constant

    # Check original and neighboring periodic images
    for i in -1:1
        x_shifted = x + i * a
        if seg.start <= x_shifted <= seg.stop
            return true
        end
    end

    return false
end

# ============================================================================
# 3D Discretization
# ============================================================================

"""
    discretize(geo::Geometry{Dim3}, resolution, property, [method])

Discretize a 3D geometry onto a grid.

# Arguments
- `geo`: The 3D geometry to discretize
- `resolution`: Tuple (Nx, Ny, Nz) of grid points
- `property`: Property to extract (:ε, :μ, :ρ, :C11, :C12, :C44, etc.)
- `method`: Discretization method (default: SimpleGrid())

# Returns
A 3D array of property values on the grid.
"""
function discretize(geo::Geometry{Dim3}, resolution::Tuple{Int,Int,Int},
                    property::Symbol, method::DiscretizationMethod=SimpleGrid())
    discretize_impl(geo, resolution, property, method)
end

function discretize_impl(geo::Geometry{Dim3}, resolution::Tuple{Int,Int,Int},
                         property::Symbol, ::SimpleGrid)
    Nx, Ny, Nz = resolution
    a1, a2, a3 = geo.lattice.vectors

    result = zeros(Float64, Nx, Ny, Nz)

    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        # Fractional coordinates [0, 1)
        fx = (i - 1) / Nx
        fy = (j - 1) / Ny
        fz = (k - 1) / Nz
        # Real space coordinates
        point = fx * a1 + fy * a2 + fz * a3
        # Get material at this point
        mat = get_material(geo, point)
        # Extract requested property
        result[i, j, k] = get_property(mat, property)
    end

    return result
end

function discretize_impl(geo::Geometry{Dim3}, resolution::Tuple{Int,Int,Int},
                         property::Symbol, method::SubpixelAverage)
    Nx, Ny, Nz = resolution
    a1, a2, a3 = geo.lattice.vectors
    ns = method.samples_per_dim

    result = zeros(Float64, Nx, Ny, Nz)

    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        # Average over subpixel samples
        val = 0.0
        for sk in 1:ns, sj in 1:ns, si in 1:ns
            fx = (i - 1 + (si - 0.5) / ns) / Nx
            fy = (j - 1 + (sj - 0.5) / ns) / Ny
            fz = (k - 1 + (sk - 0.5) / ns) / Nz
            point = fx * a1 + fy * a2 + fz * a3
            mat = get_material(geo, point)
            val += get_property(mat, property)
        end
        result[i, j, k] = val / (ns * ns * ns)
    end

    return result
end

# ============================================================================
# Convenience functions for common geometries
# ============================================================================

"""
    discretize_dielectric(geo::Geometry{Dim2}, resolution; method=SimpleGrid())

Discretize ε and μ for photonic calculations.
Returns (ε, μ) arrays.
"""
function discretize_dielectric(geo::Geometry{Dim2}, resolution;
                               method::DiscretizationMethod=SimpleGrid())
    ε = discretize(geo, resolution, :ε, method)
    μ = discretize(geo, resolution, :μ, method)
    return ε, μ
end

"""
    discretize_elastic(geo::Geometry{Dim2}, resolution; method=SimpleGrid())

Discretize elastic properties for phononic calculations.
Returns (ρ, C11, C12, C44) arrays.
"""
function discretize_elastic(geo::Geometry{Dim2}, resolution;
                            method::DiscretizationMethod=SimpleGrid())
    ρ = discretize(geo, resolution, :ρ, method)
    C11 = discretize(geo, resolution, :C11, method)
    C12 = discretize(geo, resolution, :C12, method)
    C44 = discretize(geo, resolution, :C44, method)
    return ρ, C11, C12, C44
end

# ============================================================================
# Error fallback for invalid argument types
# ============================================================================

"""
    discretize(geometry, resolution, ...)

Discretize geometry on a grid.

See concrete method signatures for detailed documentation and keyword arguments.
"""
function discretize(geometry::Any, args...; kwargs...)
    error("discretize: expected geometry::Geometry as first argument, " *
          "got $(typeof(geometry))")
end
