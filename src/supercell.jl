# Last-Modified: 2025-12-16T12:30:00+09:00

#=
Supercell construction for defect mode calculations.
=#

# ============================================================================
# Helper functions
# ============================================================================

"""
    _supercell_lattice(lat::Lattice{D}, size::NTuple{D,Int}) -> Lattice{D}

Create a supercell lattice by scaling the original lattice vectors.
"""
function _supercell_lattice(lat::Lattice{Dim1}, size::Tuple{Int})
    Nx = size[1]
    a1 = lat.vectors[1][1]  # Extract scalar from Vec1
    Lattice(Nx * a1)  # Lattice(::Real) constructor for 1D
end

function _supercell_lattice(lat::Lattice{Dim2}, size::NTuple{2,Int})
    Nx, Ny = size
    a1, a2 = lat.vectors
    Lattice(Nx * a1, Ny * a2)
end

function _supercell_lattice(lat::Lattice{Dim3}, size::NTuple{3,Int})
    Nx, Ny, Nz = size
    a1, a2, a3 = lat.vectors
    Lattice(Nx * a1, Ny * a2, Nz * a3)
end

"""
    _replicate_inclusions(inclusions, lat, size; skip_positions=[]) -> Vector{Inclusion}

Replicate inclusions over the supercell, optionally skipping certain positions.
"""
function _replicate_inclusions(
    inclusions::Vector{<:Tuple{Shape{Dim1},Material}},
    lat::Lattice{Dim1},
    size::Tuple{Int};
    skip_positions::Vector{Tuple{Int}}=Tuple{Int}[],
)
    Nx = size[1]
    a1 = lat.vectors[1]
    result = Tuple{Shape{Dim1},Material}[]

    for (shape, mat) in inclusions
        for i in 0:(Nx - 1)
            if (i,) in skip_positions
                continue
            end
            offset = i * a1[1]  # 1D: scalar offset
            push!(result, (translate(shape, offset), mat))
        end
    end

    return result
end

function _replicate_inclusions(
    inclusions::Vector{<:Tuple{Shape{Dim2},Material}},
    lat::Lattice{Dim2},
    size::NTuple{2,Int};
    skip_positions::Vector{NTuple{2,Int}}=NTuple{2,Int}[],
)
    Nx, Ny = size
    a1, a2 = lat.vectors
    result = Tuple{Shape{Dim2},Material}[]

    for (shape, mat) in inclusions
        for i in 0:(Nx - 1), j in 0:(Ny - 1)
            if (i, j) in skip_positions
                continue
            end
            offset = i * a1 + j * a2
            push!(result, (translate(shape, offset), mat))
        end
    end

    return result
end

function _replicate_inclusions(
    inclusions::Vector{<:Tuple{Shape{Dim3},Material}},
    lat::Lattice{Dim3},
    size::NTuple{3,Int};
    skip_positions::Vector{NTuple{3,Int}}=NTuple{3,Int}[],
)
    Nx, Ny, Nz = size
    a1, a2, a3 = lat.vectors
    result = Tuple{Shape{Dim3},Material}[]

    for (shape, mat) in inclusions
        for i in 0:(Nx - 1), j in 0:(Ny - 1), k in 0:(Nz - 1)
            if (i, j, k) in skip_positions
                continue
            end
            offset = i * a1 + j * a2 + k * a3
            push!(result, (translate(shape, offset), mat))
        end
    end

    return result
end

# ============================================================================
# Main API
# ============================================================================

"""
    create_supercell(geo::Geometry{D}, size::NTuple{D,Int}) -> Geometry{D}
    create_supercell(geo::Geometry{D}, size::NTuple{D,Int}; point_defects=[], modified=[])

Create a supercell by replicating the unit cell geometry.

# Arguments
- `geo`: Original unit cell geometry
- `size`: Number of replications in each direction (e.g., `(5, 5)` for 5×5 supercell)

# Keyword Arguments
- `point_defects`: Vector of cell indices to skip (create vacancies)
  - 1D: `[(i,), ...]` where `i ∈ 0:Nx-1`
  - 2D: `[(i, j), ...]` where `i ∈ 0:Nx-1`, `j ∈ 0:Ny-1`
  - 3D: `[(i, j, k), ...]` where `i ∈ 0:Nx-1`, `j ∈ 0:Ny-1`, `k ∈ 0:Nz-1`

# Returns
- New `Geometry` with enlarged lattice and replicated inclusions

# Examples
```julia
# Create a 5×5 supercell
lat = square_lattice(1.0)
air = Dielectric(1.0)
rod = Dielectric(8.9)
geo = Geometry(lat, air, [(Circle([0.5, 0.5], 0.2), rod)])

geo_super = create_supercell(geo, (5, 5))

# Create supercell with point defect at center
geo_defect = create_supercell(geo, (5, 5); point_defects=[(2, 2)])
```

See also: [`translate`](@ref), [`line_defect_positions`](@ref)
"""
function create_supercell(
    geo::Geometry{Dim1}, size::Tuple{Int}; point_defects::Vector{Tuple{Int}}=Tuple{Int}[]
)
    new_lat = _supercell_lattice(geo.lattice, size)
    new_inclusions = _replicate_inclusions(
        geo.inclusions, geo.lattice, size; skip_positions=point_defects
    )
    Geometry(new_lat, geo.background, new_inclusions)
end

function create_supercell(
    geo::Geometry{Dim2},
    size::NTuple{2,Int};
    point_defects::Vector{NTuple{2,Int}}=NTuple{2,Int}[],
)
    new_lat = _supercell_lattice(geo.lattice, size)
    new_inclusions = _replicate_inclusions(
        geo.inclusions, geo.lattice, size; skip_positions=point_defects
    )
    Geometry(new_lat, geo.background, new_inclusions)
end

function create_supercell(
    geo::Geometry{Dim3},
    size::NTuple{3,Int};
    point_defects::Vector{NTuple{3,Int}}=NTuple{3,Int}[],
)
    new_lat = _supercell_lattice(geo.lattice, size)
    new_inclusions = _replicate_inclusions(
        geo.inclusions, geo.lattice, size; skip_positions=point_defects
    )
    Geometry(new_lat, geo.background, new_inclusions)
end

# Convenience: allow integer for 1D
function create_supercell(geo::Geometry{Dim1}, n::Int; kwargs...)
    create_supercell(geo, (n,); kwargs...)
end

# ============================================================================
# Line defect helper
# ============================================================================

"""
    line_defect_positions(direction::Symbol, index::Int, size) -> Vector{NTuple}

Generate cell positions for a line defect (waveguide).

# Arguments
- `direction`: Direction of the line defect
  - 2D: `:x` (horizontal) or `:y` (vertical)
  - 3D: `:x`, `:y`, or `:z`
- `index`: Row/column index for the line defect (0-based)
- `size`: Supercell size tuple

# Returns
- Vector of cell positions to skip

# Examples
```julia
# Horizontal waveguide at row j=2 in 7×5 supercell
positions = line_defect_positions(:x, 2, (7, 5))
# Returns [(0,2), (1,2), (2,2), (3,2), (4,2), (5,2), (6,2)]

# Use with create_supercell
geo_waveguide = create_supercell(geo, (7, 5); point_defects=positions)
```
"""
function line_defect_positions(direction::Symbol, index::Int, size::NTuple{2,Int})
    Nx, Ny = size
    if direction == :x
        # Horizontal line: all cells at row j=index
        return [(i, index) for i in 0:(Nx - 1)]
    elseif direction == :y
        # Vertical line: all cells at column i=index
        return [(index, j) for j in 0:(Ny - 1)]
    else
        error("Unknown direction: $direction. Use :x or :y for 2D.")
    end
end

function line_defect_positions(direction::Symbol, index::Int, size::NTuple{3,Int})
    Nx, Ny, Nz = size
    if direction == :x
        # Line along x-axis at (y=index[1], z=index[2]) - but index is single Int
        # For simplicity, assume index is for the y-direction, z=Nz÷2
        z_idx = Nz ÷ 2
        return [(i, index, z_idx) for i in 0:(Nx - 1)]
    elseif direction == :y
        x_idx = Nx ÷ 2
        return [(x_idx, j, index) for j in 0:(Ny - 1)]
    elseif direction == :z
        x_idx = Nx ÷ 2
        return [(x_idx, index, k) for k in 0:(Nz - 1)]
    else
        error("Unknown direction: $direction. Use :x, :y, or :z for 3D.")
    end
end

"""
    line_defect_positions(direction::Symbol, indices::NTuple{2,Int}, size::NTuple{3,Int})

3D line defect with explicit (index1, index2) specification.

# Examples
```julia
# Line along x-axis at (y=2, z=3) in 7×5×5 supercell
positions = line_defect_positions(:x, (2, 3), (7, 5, 5))
```
"""
function line_defect_positions(
    direction::Symbol, indices::NTuple{2,Int}, size::NTuple{3,Int}
)
    Nx, Ny, Nz = size
    idx1, idx2 = indices
    if direction == :x
        return [(i, idx1, idx2) for i in 0:(Nx - 1)]
    elseif direction == :y
        return [(idx1, j, idx2) for j in 0:(Ny - 1)]
    elseif direction == :z
        return [(idx1, idx2, k) for k in 0:(Nz - 1)]
    else
        error("Unknown direction: $direction. Use :x, :y, or :z for 3D.")
    end
end

# Error fallback for invalid create_supercell arguments
function create_supercell(args...)
    arg_types = isempty(args) ? "()" : "(" * join(typeof.(args), ", ") * ")"
    error(
        "No matching create_supercell method for arguments: $arg_types. " *
        "Use create_supercell(geo::Geometry, size::NTuple{D,Int}; point_defects=[]).",
    )
end
