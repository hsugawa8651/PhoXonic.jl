# src/field.jl
# Field reconstruction and visualization core functions for PhoXonic.jl

#=
This file provides functions to convert eigenvectors (plane wave coefficients)
to real-space field distributions.

Key functions:
- reconstruct_field: Convert eigenvector to real-space field
- get_epsilon_field: Get permittivity distribution
- get_material_field: Get arbitrary material property distribution
- fix_phase: Normalize field phase for visualization
- field_energy: Compute energy density

Dependencies:
- Uses fourier_to_grid! from matrixfree.jl (must be included after matrixfree.jl)
=#

# =============================================================================
# reconstruct_field - Main API
# =============================================================================

"""
    reconstruct_field(solver, eigenvector; grid=nothing)

Convert an eigenvector (plane wave coefficients) to a real-space field.

# Arguments
- `solver::Solver`: The solver object containing basis and resolution info
- `eigenvector::AbstractVector`: Eigenvector from `solve_at_k_with_vectors`

# Keyword Arguments
- `grid::Union{Nothing, NTuple{D,Int}}`: Output grid size. If `nothing`, uses `solver.resolution`

# Returns
- For scalar fields (ncomponents=1): Array of dimension D
- For vector fields (ncomponents>1): Tuple of arrays

# Examples
```julia
freqs, vecs = solve_at_k_with_vectors(solver, k)
field = reconstruct_field(solver, vecs[:, 1])

# Custom grid size
field_hires = reconstruct_field(solver, vecs[:, 1]; grid=(128, 128))
```
"""
function reconstruct_field(
    solver::Solver{D,W},
    eigenvector::AbstractVector{<:Complex};
    grid::Union{Nothing,NTuple{N,Int} where N}=nothing
) where {D,W}
    resolution = isnothing(grid) ? solver.resolution : grid
    ncomp = ncomponents(solver.wave)

    if ncomp == 1
        return _reconstruct_scalar_field(solver, eigenvector, resolution)
    else
        return _reconstruct_vector_field(solver, eigenvector, resolution, ncomp)
    end
end

# -----------------------------------------------------------------------------
# 1D scalar field reconstruction
# -----------------------------------------------------------------------------
function _reconstruct_scalar_field(
    solver::Solver{Dim1,W},
    eigenvector::AbstractVector{<:Complex},
    resolution::NTuple{1,Int}
) where {W}
    N = resolution[1]
    grid = zeros(ComplexF64, N)
    fourier_to_grid!(grid, ComplexF64.(eigenvector), solver.basis, resolution)
    return grid
end

# -----------------------------------------------------------------------------
# 2D scalar field reconstruction
# -----------------------------------------------------------------------------
function _reconstruct_scalar_field(
    solver::Solver{Dim2,W},
    eigenvector::AbstractVector{<:Complex},
    resolution::NTuple{2,Int}
) where {W}
    Nx, Ny = resolution
    grid = zeros(ComplexF64, Nx, Ny)
    fourier_to_grid!(grid, ComplexF64.(eigenvector), solver.basis, resolution)
    return grid
end

# -----------------------------------------------------------------------------
# 3D scalar field reconstruction
# -----------------------------------------------------------------------------
function _reconstruct_scalar_field(
    solver::Solver{Dim3,W},
    eigenvector::AbstractVector{<:Complex},
    resolution::NTuple{3,Int}
) where {W}
    Nx, Ny, Nz = resolution
    grid = zeros(ComplexF64, Nx, Ny, Nz)
    fourier_to_grid!(grid, ComplexF64.(eigenvector), solver.basis, resolution)
    return grid
end

# -----------------------------------------------------------------------------
# 2D vector field reconstruction (ncomponents > 1)
# -----------------------------------------------------------------------------
function _reconstruct_vector_field(
    solver::Solver{Dim2,W},
    eigenvector::AbstractVector{<:Complex},
    resolution::NTuple{2,Int},
    ncomp::Int
) where {W}
    npw = solver.basis.num_pw
    @assert length(eigenvector) == ncomp * npw "Eigenvector length mismatch: expected $(ncomp * npw), got $(length(eigenvector))"

    Nx, Ny = resolution
    fields = ntuple(ncomp) do i
        start_idx = (i - 1) * npw + 1
        end_idx = i * npw
        coeffs = ComplexF64.(eigenvector[start_idx:end_idx])
        grid = zeros(ComplexF64, Nx, Ny)
        fourier_to_grid!(grid, coeffs, solver.basis, resolution)
        grid
    end

    return fields
end

# -----------------------------------------------------------------------------
# 3D vector field reconstruction (ncomponents > 1)
# -----------------------------------------------------------------------------
function _reconstruct_vector_field(
    solver::Solver{Dim3,W},
    eigenvector::AbstractVector{<:Complex},
    resolution::NTuple{3,Int},
    ncomp::Int
) where {W}
    npw = solver.basis.num_pw
    @assert length(eigenvector) == ncomp * npw "Eigenvector length mismatch"

    Nx, Ny, Nz = resolution
    fields = ntuple(ncomp) do i
        start_idx = (i - 1) * npw + 1
        end_idx = i * npw
        coeffs = ComplexF64.(eigenvector[start_idx:end_idx])
        grid = zeros(ComplexF64, Nx, Ny, Nz)
        fourier_to_grid!(grid, coeffs, solver.basis, resolution)
        grid
    end

    return fields
end

# =============================================================================
# get_epsilon_field / get_material_field
# =============================================================================

"""
    get_epsilon_field(solver)

Get the permittivity distribution on the real-space grid.

# Returns
- 1D: `Vector{Float64}`
- 2D: `Matrix{Float64}`
- 3D: `Array{Float64,3}`

# Examples
```julia
solver = Solver(TMWave(), geo, (64, 64); cutoff=7)
eps = get_epsilon_field(solver)
heatmap(eps)
```
"""
function get_epsilon_field(solver::Solver)
    mats = solver.material_arrays
    if haskey(mats, :ε)
        return copy(mats.ε)
    elseif haskey(mats, :ε_inv)
        # ε_inv only (TransverseEM, etc.)
        return 1.0 ./ mats.ε_inv
    else
        error("Solver does not have epsilon data. This may be a phononic solver.")
    end
end

"""
    get_material_field(solver, property::Symbol)

Get a material property distribution on the real-space grid.

# Arguments
- `property::Symbol`: One of `:ε`, `:μ`, `:ρ`, `:C11`, `:C12`, `:C44`

# Examples
```julia
# Phononic: get density
rho = get_material_field(solver, :ρ)

# Photonic: get permeability
mu = get_material_field(solver, :μ)
```
"""
function get_material_field(solver::Solver, property::Symbol)
    mats = solver.material_arrays
    if haskey(mats, property)
        return copy(getproperty(mats, property))
    else
        available = keys(mats)
        error("Property :$property not available. Available: $available")
    end
end

# =============================================================================
# fix_phase
# =============================================================================

"""
    fix_phase(field; method=:max)

Normalize the phase of a complex field for visualization.

After calling this function, plotting `real(field)` gives meaningful results.

# Keyword Arguments
- `method::Symbol`:
  - `:max` - Make the maximum amplitude point real and positive (default)
  - `:center` - Make the center point real and positive
  - `:mean` - Rotate to minimize imaginary part (mean phase = 0)

# Returns
- Phase-normalized field (same type as input)

# Examples
```julia
field = reconstruct_field(solver, vecs[:, 1])
field_fixed = fix_phase(field)
heatmap(real.(field_fixed))
```
"""
function fix_phase(field::AbstractArray{<:Complex}; method::Symbol=:max)
    if method == :max
        # Find maximum amplitude point
        idx = argmax(abs.(field))
        phase = angle(field[idx])
    elseif method == :center
        # Center point
        center_idx = div.(size(field) .+ 1, 2)
        phase = angle(field[center_idx...])
    elseif method == :mean
        # Mean phase (without Statistics.mean dependency)
        phase = angle(sum(field) / length(field))
    else
        error("Unknown method: $method. Use :max, :center, or :mean")
    end

    return field .* exp(-im * phase)
end

# =============================================================================
# field_energy
# =============================================================================

"""
    field_energy(solver, eigenvector)

Compute the energy density of a mode.

# Returns
- Energy density on the real-space grid (always real, non-negative)

For photonic modes: proportional to |H|² (or |E|² for TM)
For phononic modes: proportional to ρ|u|²

# Examples
```julia
energy = field_energy(solver, vecs[:, 1])
heatmap(energy)
```
"""
function field_energy(solver::Solver, eigenvector::AbstractVector{<:Complex})
    field = reconstruct_field(solver, eigenvector)
    return _compute_energy(field)
end

function _compute_energy(field::AbstractArray{<:Complex})
    return abs2.(field)
end

function _compute_energy(fields::Tuple)
    # Sum of |component|² for vector fields
    return sum(abs2.(f) for f in fields)
end

# =============================================================================
# Plot function stubs (implemented in Extensions)
# =============================================================================

"""
    plot_field(solver, eigenvector; kwargs...)

Visualize a field. Requires Plots.jl or Makie.jl.

See `PhoXonicPlotsExt` or `PhoXonicMakieExt` for implementation.

# Keyword Arguments
- `component::Symbol = :auto` - Component to plot (for vector fields)
- `quantity::Symbol = :real` - `:real`, `:imag`, `:abs`, `:phase`
- `colormap::Symbol = :viridis` - Colormap
- `title::String = "Field"` - Plot title

For 3D fields:
- `plane::Symbol = :xy` - Slice plane (`:xy`, `:yz`, `:xz`)
- `slice::Real = 0.5` - Slice position (0 to 1)
"""
function plot_field end

"""
    plot_epsilon(solver; kwargs...)

Visualize permittivity distribution. Requires Plots.jl or Makie.jl.
"""
function plot_epsilon end

"""
    plot_field!(fig, solver, eigenvector; kwargs...)

Add field to existing figure.
"""
function plot_field! end
