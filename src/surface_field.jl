"""
    surface_field(solver, eigenvector; kwargs...)

Plot a 3D surface representation of a 2D field as a Makie `Figure`.

Implementations live in `PhoXonicMakieExt` and are loaded when `Makie` is
available in the active environment.

# Arguments
- `solver::Solver{Dim2}`: the 2D solver used to compute eigenvectors.
- `eigenvector::AbstractVector`: the eigenvector to visualize.

# Keyword arguments
See `ext/PhoXonicMakieExt.jl` for the full keyword list (`component`,
`quantity`, `colormap`, `title`, `size`, etc.).
"""
function surface_field end
