# Band structure calculation using Transfer Matrix Method
#
# For a periodic structure, Bloch's theorem gives:
#   Tr(M) = 2*cos(k*a)
# where M is the transfer matrix for one period (unit cell) and a is the period.
#
# For |Tr(M)| <= 2: propagating wave (real k)
# For |Tr(M)| > 2:  evanescent wave (bandgap, complex k)
#
# To find ω(k), we sweep ω and find where Tr(M(ω)) = 2*cos(k*a).

using LinearAlgebra
using StaticArrays

"""
    bloch_k(M::Matrix, a::Real) -> ComplexF64

Compute the Bloch wavenumber k from a transfer matrix M.

For |Tr(M)| <= 2, k is real (propagating).
For |Tr(M)| > 2, k is complex (bandgap).

# Arguments
- `M`: 2×2 transfer matrix for one period
- `a`: period length

# Returns
- `k`: Bloch wavenumber (may be complex in bandgap)
"""
function bloch_k(M::AbstractMatrix, a::Real)
    trace = tr(M)
    # k = (1/a) * acos(Tr(M)/2)
    # acos returns complex for |arg| > 1
    return acos(trace / 2) / a
end

"""
    tmm_bandstructure(solver::TMMSolver; k_points=51, bands=1:5, ω_max=nothing)

Compute photonic band structure using the Transfer Matrix Method.

Uses Bloch's theorem: for a periodic structure, the dispersion relation
ω(k) is found by solving Tr(M(ω)) = 2*cos(k*a).

# Arguments
- `solver::TMMSolver`: TMM solver with unit cell definition
- `k_points::Int`: number of k-points from 0 to π/a (default: 51)
- `bands`: range of bands to compute (default: 1:5)
- `ω_max`: maximum frequency to search (default: auto)

# Returns
- `BandStructure`: band structure result compatible with PWE output

# Example
```julia
unit_cell = [Layer(Dielectric(9.0), 0.25), Layer(Dielectric(1.0), 0.75)]
ml = Multilayer(unit_cell, Dielectric(1.0), Dielectric(1.0))
solver = TMMSolver(Photonic1D(), ml)

bands = tmm_bandstructure(solver; k_points=51, bands=1:5)
```
"""
function tmm_bandstructure(solver::TMMSolver{Photonic1D};
                           k_points::Int=51,
                           bands::Union{Int,UnitRange{Int}}=1:5,
                           ω_max::Union{Nothing,Real}=nothing)
    ml = solver.structure

    # Period = total thickness of unit cell
    a = total_thickness(ml)

    # Convert bands to range if integer
    band_range = bands isa Int ? (1:bands) : bands
    nbands = length(band_range)

    # k-points from 0 to π/a (first Brillouin zone)
    k_values = range(0.0, π/a, length=k_points)

    # Estimate ω_max if not provided
    # For photonic crystal, ω ~ c*k/n_eff, use n_min for upper bound
    if ω_max === nothing
        n_min = minimum(_get_n(layer.material) for layer in ml.layers)
        n_min = real(n_min)
        ω_max = 2π * maximum(band_range) / (a * n_min) * 1.5
    end

    # Storage for frequencies
    frequencies = zeros(k_points, nbands)

    # For each k, find frequencies where Tr(M(ω)) = 2*cos(k*a)
    for (ik, k) in enumerate(k_values)
        target = 2 * cos(k * a)
        ω_list = find_band_frequencies(ml, target, nbands, ω_max, a)

        for (ib, ω) in enumerate(ω_list)
            if ib <= nbands
                frequencies[ik, ib] = ω
            end
        end
    end

    # Construct k-path (1D: just a line from 0 to π/a)
    kpoints = [SVector{1}(k) for k in k_values]

    # Compute distances along the path
    distances = collect(k_values)  # For 1D, distance = k value

    # Labels for high-symmetry points
    labels = [(1, "Γ"), (k_points, "X")]

    # Return BandStructure compatible with PWE
    return BandStructure{1}(kpoints, distances, frequencies, labels)
end

"""
    find_band_frequencies(ml, target, nbands, ω_max, a) -> Vector{Float64}

Find frequencies where Tr(M(ω)) = target.

Sweeps frequency from 0 to ω_max and detects crossings and touch points.
For zone boundaries (target = ±2), also detects where trace touches target.
"""
function find_band_frequencies(ml::Multilayer, target::Real,
                               nbands::Int, ω_max::Real, a::Real)
    # Number of points to sweep (higher = more accurate)
    nω = 2000
    ω_values = range(1e-6, ω_max, length=nω)

    # Compute Tr(M) for each ω
    # Note: ω = 2πc/λ, for normalized units (c=1), λ = 2π/ω
    trace_values = zeros(nω)
    for (i, ω) in enumerate(ω_values)
        λ = 2π / ω  # wavelength from angular frequency
        M = system_matrix(ml, λ)
        trace_values[i] = real(tr(M))
    end

    # Find crossings and touch points where trace equals target
    frequencies = Float64[]
    tol = 0.01  # Tolerance for touch detection

    # At k=0 (target=2), ω=0 is always a solution (first band at Γ)
    if abs(target - 2.0) < 1e-10
        push!(frequencies, 0.0)
    end

    for i in 1:(nω-1)
        t1, t2 = trace_values[i], trace_values[i+1]

        # Check if target is crossed
        if (t1 - target) * (t2 - target) < 0
            # Linear interpolation to find crossing point
            ω1, ω2 = ω_values[i], ω_values[i+1]
            ω_cross = ω1 + (target - t1) / (t2 - t1) * (ω2 - ω1)
            push!(frequencies, ω_cross)
        # Check for touch points at zone boundaries (target ≈ ±2)
        elseif abs(target) >= 2.0 - tol && i > 1 && i < nω - 1
            t0 = trace_values[i-1]
            # Detect local extremum near target
            # At zone boundary, trace touches ±2 from inside propagating region
            if target > 0  # target ≈ +2 (Γ point)
                # Look for local maximum approaching +2
                if t1 > t0 && t1 > t2 && abs(t1 - target) < tol
                    push!(frequencies, ω_values[i])
                end
            else  # target ≈ -2 (X point for some structures)
                # Look for local minimum approaching -2
                if t1 < t0 && t1 < t2 && abs(t1 - target) < tol
                    push!(frequencies, ω_values[i])
                end
            end
        end

        if length(frequencies) >= nbands
            break
        end
    end

    # Sort and remove duplicates (within tolerance)
    sort!(frequencies)
    unique_freqs = Float64[]
    for ω in frequencies
        if isempty(unique_freqs) || ω - unique_freqs[end] > ω_max / nω * 5
            push!(unique_freqs, ω)
        end
    end

    # Pad with zeros if not enough bands found
    while length(unique_freqs) < nbands
        push!(unique_freqs, 0.0)
    end

    return unique_freqs[1:nbands]
end

# ============================================================================
# Phononic (Longitudinal1D) band structure
# ============================================================================

"""
    tmm_bandstructure(solver::TMMSolver{Longitudinal1D}; k_points=51, bands=1:5, ω_max=nothing)

Compute phononic band structure using the Transfer Matrix Method.

Uses Bloch's theorem: for a periodic structure, the dispersion relation
ω(k) is found by solving Tr(M(ω)) = 2*cos(k*a).

# Arguments
- `solver::TMMSolver{Longitudinal1D}`: TMM solver with unit cell definition
- `k_points::Int`: number of k-points from 0 to π/a (default: 51)
- `bands`: range of bands to compute (default: 1:5)
- `ω_max`: maximum angular frequency to search (default: auto)

# Returns
- `BandStructure`: band structure result compatible with PWE output

# Example
```julia
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
unit_cell = [Layer(steel, 0.0025), Layer(epoxy, 0.0075)]
ml = Multilayer(unit_cell, epoxy, epoxy)
solver = TMMSolver(Longitudinal1D(), ml)

bands = tmm_bandstructure(solver; k_points=51, bands=1:5)
```
"""
function tmm_bandstructure(solver::TMMSolver{Longitudinal1D};
                           k_points::Int=51,
                           bands::Union{Int,UnitRange{Int}}=1:5,
                           ω_max::Union{Nothing,Real}=nothing)
    ml = solver.structure

    # Period = total thickness of unit cell
    a = total_thickness(ml)

    # Convert bands to range if integer
    band_range = bands isa Int ? (1:bands) : bands
    nbands = length(band_range)

    # k-points from 0 to π/a (first Brillouin zone)
    k_values = range(0.0, π/a, length=k_points)

    # Estimate ω_max if not provided
    # For phononic crystal, ω ~ c*k, use c_max for upper bound
    if ω_max === nothing
        c_max = maximum(longitudinal_velocity(layer.material) for layer in ml.layers)
        ω_max = 2π * maximum(band_range) / a * c_max * 1.5
    end

    # Storage for frequencies
    frequencies = zeros(k_points, nbands)

    # For each k, find frequencies where Tr(M(ω)) = 2*cos(k*a)
    for (ik, k) in enumerate(k_values)
        target = 2 * cos(k * a)
        ω_list = find_band_frequencies_acoustic(ml, target, nbands, ω_max, a)

        for (ib, ω) in enumerate(ω_list)
            if ib <= nbands
                frequencies[ik, ib] = ω
            end
        end
    end

    # Construct k-path (1D: just a line from 0 to π/a)
    kpoints = [SVector{1}(k) for k in k_values]

    # Compute distances along the path
    distances = collect(k_values)

    # Labels for high-symmetry points
    labels = [(1, "Γ"), (k_points, "X")]

    return BandStructure{1}(kpoints, distances, frequencies, labels)
end

"""
    find_band_frequencies_acoustic(ml, target, nbands, ω_max, a) -> Vector{Float64}

Find frequencies where Tr(M_acoustic(ω)) = target for phononic crystals.
"""
function find_band_frequencies_acoustic(ml::Multilayer{<:ElasticMaterial}, target::Real,
                                        nbands::Int, ω_max::Real, a::Real)
    # Number of points to sweep
    nω = 2000
    ω_values = range(1e-6, ω_max, length=nω)

    # Compute Tr(M) for each ω
    # Note: for phononic, ω = 2π*c/λ, so λ = 2π*c/ω
    # We need average velocity for wavelength conversion
    c_avg = sum(longitudinal_velocity(layer.material) * thickness(layer)
                for layer in ml.layers) / a

    trace_values = zeros(nω)
    for (i, ω) in enumerate(ω_values)
        λ = 2π * c_avg / ω  # wavelength from angular frequency
        M = system_matrix_acoustic(ml, λ)
        trace_values[i] = real(tr(M))
    end

    # Find crossings and touch points
    frequencies = Float64[]
    tol = 0.01

    # At k=0 (target=2), ω=0 is always a solution
    if abs(target - 2.0) < 1e-10
        push!(frequencies, 0.0)
    end

    for i in 1:(nω-1)
        t1, t2 = trace_values[i], trace_values[i+1]

        # Check if target is crossed
        if (t1 - target) * (t2 - target) < 0
            ω1, ω2 = ω_values[i], ω_values[i+1]
            ω_cross = ω1 + (target - t1) / (t2 - t1) * (ω2 - ω1)
            push!(frequencies, ω_cross)
        elseif abs(target) >= 2.0 - tol && i > 1 && i < nω - 1
            t0 = trace_values[i-1]
            if target > 0
                if t1 > t0 && t1 > t2 && abs(t1 - target) < tol
                    push!(frequencies, ω_values[i])
                end
            else
                if t1 < t0 && t1 < t2 && abs(t1 - target) < tol
                    push!(frequencies, ω_values[i])
                end
            end
        end

        if length(frequencies) >= nbands
            break
        end
    end

    # Sort and remove duplicates
    sort!(frequencies)
    unique_freqs = Float64[]
    for ω in frequencies
        if isempty(unique_freqs) || ω - unique_freqs[end] > ω_max / nω * 5
            push!(unique_freqs, ω)
        end
    end

    # Pad with zeros if not enough bands found
    while length(unique_freqs) < nbands
        push!(unique_freqs, 0.0)
    end

    return unique_freqs[1:nbands]
end
