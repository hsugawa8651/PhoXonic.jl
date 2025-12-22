# Last-Modified: 2025-12-16T12:00:00+09:00

#=
Geometric shapes for defining crystal structures.
=#

"""
Abstract type for geometric shapes in D dimensions.
"""
abstract type Shape{D<:Dimension} end

# ============================================================================
# 2D Shapes
# ============================================================================

"""
    Circle(center, radius)

A circle in 2D.

# Examples
```julia
c = Circle(Vec2(0, 0), 0.2)
c = Circle([0.0, 0.0], 0.2)
```
"""
struct Circle <: Shape{Dim2}
    center::Vec2
    radius::Float64
end

Circle(center::AbstractVector, radius::Real) = Circle(Vec2(center...), Float64(radius))

"""
    in(point, circle::Circle)

Check if a point is inside the circle.
"""
function Base.in(point::Vec2, c::Circle)
    norm(point - c.center) <= c.radius
end

function Base.in(point::AbstractVector, c::Circle)
    Vec2(point...) in c
end

"""
    Rectangle(center, size)

An axis-aligned rectangle in 2D.

# Examples
```julia
r = Rectangle(Vec2(0, 0), Vec2(0.5, 0.3))
```
"""
struct Rectangle <: Shape{Dim2}
    center::Vec2
    size::Vec2  # (width, height)
end

function Rectangle(center::AbstractVector, size::AbstractVector)
    Rectangle(Vec2(center...), Vec2(size...))
end

"""
    in(point, rect::Rectangle)

Check if a point is inside the rectangle.
"""
function Base.in(point::Vec2, r::Rectangle)
    dx = abs(point[1] - r.center[1])
    dy = abs(point[2] - r.center[2])
    dx <= r.size[1] / 2 && dy <= r.size[2] / 2
end

function Base.in(point::AbstractVector, r::Rectangle)
    Vec2(point...) in r
end

"""
    Polygon(vertices)

A polygon in 2D defined by its vertices (in order).
"""
struct Polygon <: Shape{Dim2}
    vertices::Vector{Vec2}
end

Polygon(vertices::Vector{<:AbstractVector}) = Polygon([Vec2(v...) for v in vertices])

"""
    in(point, poly::Polygon)

Check if a point is inside the polygon using ray casting algorithm.
"""
function Base.in(point::Vec2, poly::Polygon)
    n = length(poly.vertices)
    inside = false
    j = n
    for i in 1:n
        vi = poly.vertices[i]
        vj = poly.vertices[j]
        if ((vi[2] > point[2]) != (vj[2] > point[2])) &&
            (point[1] < (vj[1] - vi[1]) * (point[2] - vi[2]) / (vj[2] - vi[2]) + vi[1])
            inside = !inside
        end
        j = i
    end
    inside
end

# ============================================================================
# 3D Shapes
# ============================================================================

"""
    Sphere(center, radius)

A sphere in 3D.
"""
struct Sphere <: Shape{Dim3}
    center::Vec3
    radius::Float64
end

Sphere(center::AbstractVector, radius::Real) = Sphere(Vec3(center...), Float64(radius))

"""
    in(point, sphere::Sphere)

Check if a point is inside the sphere.
"""
function Base.in(point::Vec3, s::Sphere)
    norm(point - s.center) <= s.radius
end

function Base.in(point::AbstractVector, s::Sphere)
    Vec3(point...) in s
end

"""
    Cylinder(center, radius, height, axis)

A cylinder in 3D with specified center, radius, height, and axis direction.
Default axis is along z.
"""
struct Cylinder <: Shape{Dim3}
    center::Vec3
    radius::Float64
    height::Float64
    axis::Vec3  # Normalized axis direction
end

function Cylinder(center::Vec3, radius::Real, height::Real, axis::Vec3=Vec3(0, 0, 1))
    Cylinder(center, Float64(radius), Float64(height), normalize(axis))
end

function Cylinder(
    center::AbstractVector, radius::Real, height::Real, axis::AbstractVector=Vec3(0, 0, 1)
)
    Cylinder(Vec3(center...), Float64(radius), Float64(height), normalize(Vec3(axis...)))
end

"""
    in(point, cyl::Cylinder)

Check if a point is inside the cylinder.
"""
function Base.in(point::Vec3, cyl::Cylinder)
    # Project point onto axis
    d = point - cyl.center
    z = dot(d, cyl.axis)
    # Check height
    if abs(z) > cyl.height / 2
        return false
    end
    # Check radial distance
    r_vec = d - z * cyl.axis
    norm(r_vec) <= cyl.radius
end

"""
    Slab(z_min, z_max)

An infinite slab in the xy-plane, bounded by z_min ≤ z ≤ z_max.

# Examples
```julia
# Si layer from z=0.3 to z=0.8
slab = Slab(0.3, 0.8)
```
"""
struct Slab <: Shape{Dim3}
    z_min::Float64
    z_max::Float64
end

function Slab(z_min::Real, z_max::Real)
    @assert z_min < z_max "z_min must be less than z_max"
    Slab(Float64(z_min), Float64(z_max))
end

"""
    in(point, slab::Slab)

Check if a point is inside the slab (z_min ≤ z ≤ z_max).
"""
function Base.in(point::Vec3, slab::Slab)
    slab.z_min <= point[3] <= slab.z_max
end

function Base.in(point::AbstractVector, slab::Slab)
    Vec3(point...) in slab
end

# ============================================================================
# 1D Shapes
# ============================================================================

"""
    Segment(start, stop)

A line segment in 1D.
"""
struct Segment <: Shape{Dim1}
    start::Float64
    stop::Float64
end

"""
    in(point, seg::Segment)

Check if a point is inside the segment.
"""
function Base.in(point::Real, seg::Segment)
    seg.start <= point <= seg.stop
end

function Base.in(point::Vec1, seg::Segment)
    point[1] in seg
end

# ============================================================================
# Shape Translation (for supercell construction)
# ============================================================================

"""
    translate(shape::Shape, offset) -> Shape

Create a new shape translated by the given offset vector.

The offset type should match the shape dimension:
- 1D shapes: `offset::Real`
- 2D shapes: `offset::Vec2` or `AbstractVector` (length 2)
- 3D shapes: `offset::Vec3` or `AbstractVector` (length 3)

# Examples
```julia
# 2D
c = Circle([0.0, 0.0], 0.2)
c2 = translate(c, [1.0, 2.0])  # center at (1, 2)

# 3D
s = Sphere([0.0, 0.0, 0.0], 0.3)
s2 = translate(s, [1.0, 0.0, 0.5])  # center at (1, 0, 0.5)
```
"""
function translate end

# ----------------------------------------------------------------------------
# 2D Shape translations
# ----------------------------------------------------------------------------

translate(c::Circle, offset::Vec2) = Circle(c.center + offset, c.radius)
translate(c::Circle, offset::AbstractVector) = translate(c, Vec2(offset...))

translate(r::Rectangle, offset::Vec2) = Rectangle(r.center + offset, r.size)
translate(r::Rectangle, offset::AbstractVector) = translate(r, Vec2(offset...))

translate(p::Polygon, offset::Vec2) = Polygon([v + offset for v in p.vertices])
translate(p::Polygon, offset::AbstractVector) = translate(p, Vec2(offset...))

# ----------------------------------------------------------------------------
# 3D Shape translations
# ----------------------------------------------------------------------------

translate(s::Sphere, offset::Vec3) = Sphere(s.center + offset, s.radius)
translate(s::Sphere, offset::AbstractVector) = translate(s, Vec3(offset...))

function translate(c::Cylinder, offset::Vec3)
    Cylinder(c.center + offset, c.radius, c.height, c.axis)
end
translate(c::Cylinder, offset::AbstractVector) = translate(c, Vec3(offset...))

# Slab: xy-infinite, only z-translation is meaningful
translate(s::Slab, offset::Vec3) = Slab(s.z_min + offset[3], s.z_max + offset[3])
translate(s::Slab, offset::AbstractVector) = translate(s, Vec3(offset...))

# ----------------------------------------------------------------------------
# 1D Shape translations
# ----------------------------------------------------------------------------

translate(s::Segment, offset::Real) = Segment(s.start + offset, s.stop + offset)
