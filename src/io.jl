# Last-Modified: 2025-12-14T17:50:17+09:00

#=
I/O utilities for PhoXonic.jl

Provides functions to save and load band structures and other results using JLD2.

Usage:
    using PhoXonic

    bands = compute_bands(solver, kpath)
    save_bands("bands.jld2", bands)

    # Later...
    bands_loaded = load_bands("bands.jld2")
=#

using JLD2
using Dates: now

# Version for file format compatibility
const IO_VERSION = v"1.0.0"

# ============================================================================
# Band Structure I/O
# ============================================================================

"""
    save_bands(filename, bands::BandStructure; compress=false, metadata=Dict())

Save a band structure to a JLD2 file.

# Arguments
- `filename`: Path to the output file (should end with `.jld2`)
- `bands`: BandStructure object from `compute_bands()`

# Keyword Arguments
- `compress`: Enable compression (default: false)
- `metadata`: Additional metadata dictionary

# Example
```julia
using PhoXonic

solver = Solver(TEWave(), geo, (64, 64))
bands = compute_bands(solver, kpath)
save_bands("te_bands.jld2", bands)
```

# Saved Data
The file contains:
- `frequencies`: Matrix of frequencies (nkpoints × nbands)
- `kpoints`: Vector of k-point coordinates
- `distances`: Cumulative distances along k-path
- `labels`: High-symmetry point labels and indices
- `metadata`: Additional information (version, date, etc.)
"""
function save_bands(
    filename::AbstractString,
    bands::BandStructure;
    compress::Bool=false,
    metadata::Dict=Dict{String,Any}(),
)
    # Prepare metadata
    meta = merge(
        Dict{String,Any}(
            "created" => string(now()),
            "nkpoints" => length(bands.distances),
            "nbands" => size(bands.frequencies, 2),
        ),
        metadata,
    )

    # Save using JLD2
    jldopen(filename, "w"; compress=compress) do file
        file["version"] = string(IO_VERSION)
        file["type"] = "BandStructure"
        file["frequencies"] = bands.frequencies
        file["kpoints"] = bands.kpoints
        file["distances"] = bands.distances
        file["labels"] = bands.labels
        file["metadata"] = meta
    end

    return filename
end

"""
    load_bands(filename) -> BandStructure

Load a band structure from a JLD2 file.

# Arguments
- `filename`: Path to the input file

# Returns
A `BandStructure` object.

# Example
```julia
using PhoXonic

bands = load_bands("te_bands.jld2")
plot_bands(bands)
```
"""
function load_bands(filename::AbstractString)
    file_data = load(filename)

    # Check version compatibility
    if haskey(file_data, "version")
        file_version = VersionNumber(file_data["version"])
        if file_version > IO_VERSION
            @warn "File was created with a newer version of PhoXonic.jl" file_version IO_VERSION
        end
    end

    # Reconstruct BandStructure
    # Note: BandStructure{D}(kpoints, distances, frequencies, labels)
    kpoints = file_data["kpoints"]
    distances = file_data["distances"]
    frequencies = file_data["frequencies"]
    labels = file_data["labels"]

    return BandStructure(kpoints, distances, frequencies, labels)
end

# ============================================================================
# Mode I/O
# ============================================================================

"""
    save_modes(filename, modes; k=nothing, frequencies=nothing, compress=false, metadata=Dict())

Save eigenmodes to a JLD2 file.

# Arguments
- `filename`: Output file path
- `modes`: Eigenmode matrix (N × nbands) from `solve()`

# Keyword Arguments
- `k`: Wave vector (optional, for reference)
- `frequencies`: Corresponding frequencies (optional)
- `compress`: Enable compression (default: false, requires CodecZlib.jl)
- `metadata`: Additional metadata dictionary

# Example
```julia
using PhoXonic

ω, modes = solve(solver, k; bands=1:5)
save_modes("modes_at_k.jld2", modes; k=k, frequencies=ω)
```
"""
function save_modes(
    filename::AbstractString,
    modes::AbstractMatrix;
    k=nothing,
    frequencies=nothing,
    compress::Bool=false,
    metadata::Dict=Dict{String,Any}(),
)
    # Prepare metadata
    meta = merge(
        Dict{String,Any}("created" => string(now()), "size" => size(modes)), metadata
    )

    # Save using JLD2
    jldopen(filename, "w"; compress=compress) do file
        file["version"] = string(IO_VERSION)
        file["type"] = "Modes"
        file["modes"] = modes
        file["metadata"] = meta
        if k !== nothing
            file["k"] = k
        end
        if frequencies !== nothing
            file["frequencies"] = frequencies
        end
    end

    return filename
end

"""
    load_modes(filename) -> NamedTuple

Load eigenmodes from a JLD2 file.

# Returns
A NamedTuple with:
- `modes`: The eigenmode matrix
- `k`: Wave vector (if saved)
- `frequencies`: Frequencies (if saved)
- `metadata`: Additional metadata
"""
function load_modes(filename::AbstractString)
    file_data = load(filename)

    result = (
        modes=file_data["modes"],
        k=get(file_data, "k", nothing),
        frequencies=get(file_data, "frequencies", nothing),
        metadata=get(file_data, "metadata", Dict()),
    )

    return result
end

# ============================================================================
# Material Arrays I/O
# ============================================================================

# Multiple dispatch helpers for saving material arrays based on wave type
# Uses abstract type fallback + concrete type override for 3D

# Error fallback for unsupported wave types
function _save_material_arrays!(file, wave, mats)
    throw(ArgumentError("Unsupported wave type for save_epsilon: $(typeof(wave))"))
end

# Fallback for all PhotonicWave (TEWave, TMWave, Photonic1D)
function _save_material_arrays!(file, ::PhotonicWave, mats)
    haskey(mats, :ε) && (file["epsilon"] = mats.ε)
    haskey(mats, :ε_inv) && (file["epsilon_inv"] = mats.ε_inv)
    haskey(mats, :μ) && (file["mu"] = mats.μ)
    haskey(mats, :μ_inv) && (file["mu_inv"] = mats.μ_inv)
    return nothing
end

# Override for FullVectorEM (3D Photonic - different field subset)
function _save_material_arrays!(file, ::FullVectorEM, mats)
    haskey(mats, :ε_inv) && (file["epsilon_inv"] = mats.ε_inv)
    haskey(mats, :μ) && (file["mu"] = mats.μ)
    return nothing
end

# Fallback for all PhononicWave (SHWave, PSVWave, Longitudinal1D)
function _save_material_arrays!(file, ::PhononicWave, mats)
    haskey(mats, :C44) && (file["C44"] = mats.C44)
    haskey(mats, :C11) && (file["C11"] = mats.C11)
    haskey(mats, :C12) && (file["C12"] = mats.C12)
    haskey(mats, :ρ) && (file["rho"] = mats.ρ)
    return nothing
end

# Override for FullElastic (3D Phononic - unconditional save)
function _save_material_arrays!(file, ::FullElastic, mats)
    file["C11"] = mats.C11
    file["C12"] = mats.C12
    file["C44"] = mats.C44
    file["rho"] = mats.ρ
    return nothing
end

"""
    save_epsilon(filename, solver::Solver; compress=false, metadata=Dict())

Save discretized material arrays to a JLD2 file.

For photonic crystals, saves:
- `epsilon`: Dielectric constant
- `epsilon_inv`: Inverse dielectric constant
- `mu`: Permeability

For phononic crystals, saves:
- `C11`, `C12`, `C44`: Elastic constants (Voigt notation)
- `rho`: Density

# Example
```julia
using PhoXonic

solver = Solver(TEWave(), geo, (64, 64))
save_epsilon("epsilon.jld2", solver)

# For phononic
solver_ph = Solver(FullElastic(), geo_phononic, (16, 16, 16))
save_epsilon("elastic_constants.jld2", solver_ph)
```
"""
function save_epsilon(
    filename::AbstractString,
    solver::Solver;
    compress::Bool=false,
    metadata::Dict=Dict{String,Any}(),
)
    mats = solver.material_arrays
    wave = solver.wave

    # Prepare metadata
    meta = merge(
        Dict{String,Any}("created" => string(now()), "wave_type" => string(typeof(wave))),
        metadata,
    )

    # Save using JLD2
    jldopen(filename, "w"; compress=compress) do file
        file["version"] = string(IO_VERSION)
        file["type"] = "MaterialArrays"
        file["resolution"] = solver.resolution
        file["metadata"] = meta

        # Save based on wave type (uses multiple dispatch)
        _save_material_arrays!(file, wave, mats)
    end

    return filename
end

"""
    load_epsilon(filename) -> NamedTuple

Load discretized material arrays from a JLD2 file.

# Returns
A NamedTuple with the material arrays and metadata.

For photonic data: `:ε`, `:ε_inv`, `:μ`, `:μ_inv`
For phononic data: `:C11`, `:C12`, `:C44`, `:ρ`
Plus: `:resolution`, `:wave_type`, `:metadata`
"""
function load_epsilon(filename::AbstractString)
    file_data = load(filename)

    # Build result NamedTuple dynamically
    result = Dict{Symbol,Any}()

    # Photonic fields
    for (file_key, julia_key) in
        [("epsilon", :ε), ("epsilon_inv", :ε_inv), ("mu", :μ), ("mu_inv", :μ_inv)]
        if haskey(file_data, file_key)
            result[julia_key] = file_data[file_key]
        end
    end

    # Phononic fields
    for (file_key, julia_key) in [("C11", :C11), ("C12", :C12), ("C44", :C44), ("rho", :ρ)]
        if haskey(file_data, file_key)
            result[julia_key] = file_data[file_key]
        end
    end

    # Metadata
    result[:resolution] = get(file_data, "resolution", nothing)
    result[:wave_type] = get(get(file_data, "metadata", Dict()), "wave_type", nothing)
    result[:metadata] = get(file_data, "metadata", Dict())

    return NamedTuple(result)
end
