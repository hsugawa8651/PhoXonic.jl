# TMM Solver for transmission/reflection spectrum calculations

"""
    TMMSolver{W<:WaveType, M<:Material}

Solver for Transfer Matrix Method calculations.

# Fields
- `wavetype::W`: Wave type (Photonic1D, Longitudinal1D)
- `structure::Multilayer{M}`: The multilayer structure

# Example
```julia
air = Dielectric(1.0)
glass = Dielectric(2.25)
layer = Layer(glass, 100.0)
ml = Multilayer([layer], air, air)

solver = TMMSolver(Photonic1D(), ml)
result = tmm_spectrum(solver, 600.0)
```
"""
struct TMMSolver{W<:WaveType,M<:Material}
    wavetype::W
    structure::Multilayer{M}
end

"""
    TMMResult

Result of TMM spectrum calculation.

# Fields
- `R::Float64`: Power reflectance (0 to 1)
- `T::Float64`: Power transmittance (0 to 1)
- `A::Float64`: Absorption (0 to 1), only nonzero for lossy materials
- `r::ComplexF64`: Amplitude reflection coefficient
- `t::ComplexF64`: Amplitude transmission coefficient
"""
struct TMMResult
    R::Float64
    T::Float64
    A::Float64
    r::ComplexF64
    t::ComplexF64
end

"""
    tmm_spectrum(solver::TMMSolver, λ::Real; angle=0.0, polarization=:TE) -> TMMResult

Calculate transmission and reflection at a single wavelength.

# Arguments
- `solver::TMMSolver`: The TMM solver
- `λ::Real`: Wavelength (same units as layer thicknesses)
- `angle::Real`: Incident angle in radians (default: 0.0 for normal incidence)
- `polarization::Symbol`: `:TE` (s-wave) or `:TM` (p-wave) for oblique incidence

# Returns
- `TMMResult`: Contains R, T, A (absorption), r, t

# Example
```julia
# Normal incidence
result = tmm_spectrum(solver, 600.0)

# Oblique incidence at 45°
result = tmm_spectrum(solver, 600.0; angle=π/4, polarization=:TM)
```
"""
function tmm_spectrum(
    solver::TMMSolver{Photonic1D}, λ::Real; angle::Real=0.0, polarization::Symbol=:TE
)
    ml = solver.structure

    # Get refractive indices
    n_inc = _get_n(ml.incident)
    n_sub = _get_n(ml.substrate)

    # Get system matrix (normal or oblique)
    if angle ≈ 0.0
        M = system_matrix(ml, λ)
    else
        M = system_matrix_oblique(ml, λ, angle, polarization)
    end

    # Extract transmission and reflection coefficients
    # M * [1; r] = [t; 0] in reverse direction
    # [E_inc^+; E_inc^-] = M * [E_sub^+; 0]
    # So: 1 = M[1,1]*t, r = M[2,1]*t
    # Therefore: t = 1/M[1,1], r = M[2,1]/M[1,1]
    t = 1 / M[1, 1]
    r = M[2, 1] / M[1, 1]

    # Power coefficients
    R = abs2(r)

    # Transmission needs refractive index and angle correction
    # For oblique incidence: T = (n_sub cos θ_sub / n_inc cos θ_inc) |t|²
    if angle ≈ 0.0
        n_inc_r = real(n_inc)
        n_sub_r = real(n_sub)
        T = (n_sub_r / n_inc_r) * abs2(t)
    else
        θ_sub = snell_angle(n_inc, n_sub, angle)
        # Use real parts for energy flow calculation
        cos_inc = cos(angle)
        cos_sub = real(cos(θ_sub))
        n_inc_r = real(n_inc)
        n_sub_r = real(n_sub)
        T = (n_sub_r * cos_sub) / (n_inc_r * cos_inc) * abs2(t)
    end

    # Absorption = 1 - R - T
    A = max(0.0, 1.0 - R - T)

    return TMMResult(R, T, A, r, t)
end

"""
    tmm_spectrum(solver::TMMSolver, λ_values::AbstractVector; angle=0.0, polarization=:TE) -> (R, T)

Calculate transmission and reflection spectrum over multiple wavelengths.

# Arguments
- `solver::TMMSolver`: The TMM solver
- `λ_values::AbstractVector`: Array of wavelengths
- `angle::Real`: Incident angle in radians (default: 0.0)
- `polarization::Symbol`: `:TE` or `:TM` (default: `:TE`)

# Returns
- `R::Vector{Float64}`: Reflectance spectrum
- `T::Vector{Float64}`: Transmittance spectrum

# Example
```julia
λ_values = 400:10:800
R, T = tmm_spectrum(solver, λ_values)

# Oblique incidence
R_te, T_te = tmm_spectrum(solver, λ_values; angle=π/6, polarization=:TE)
```
"""
function tmm_spectrum(
    solver::TMMSolver, λ_values::AbstractVector; angle::Real=0.0, polarization::Symbol=:TE
)
    n = length(λ_values)
    R = Vector{Float64}(undef, n)
    T = Vector{Float64}(undef, n)

    for (i, λ) in enumerate(λ_values)
        result = tmm_spectrum(solver, λ; angle=angle, polarization=polarization)
        R[i] = result.R
        T[i] = result.T
    end

    return R, T
end

# ============================================================================
# Phononic (Longitudinal1D) TMM solver
# ============================================================================

"""
    tmm_spectrum(solver::TMMSolver{Longitudinal1D}, λ::Real) -> TMMResult

Calculate transmission and reflection for longitudinal elastic waves.

# Arguments
- `solver::TMMSolver{Longitudinal1D}`: The TMM solver for elastic waves
- `λ::Real`: Wavelength (same units as layer thicknesses)

# Returns
- `TMMResult`: Contains R, T, A (=0 for lossless), r, t

# Example
```julia
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
layer = Layer(steel, 0.001)
ml = Multilayer([layer], epoxy, epoxy)

solver = TMMSolver(Longitudinal1D(), ml)
result = tmm_spectrum(solver, 0.01)  # λ = 10 mm
```
"""
function tmm_spectrum(solver::TMMSolver{Longitudinal1D}, λ::Real)
    ml = solver.structure

    # Get system matrix
    M = system_matrix_acoustic(ml, λ)

    # Get impedances for energy calculation
    Z_inc = acoustic_impedance(ml.incident)
    Z_sub = acoustic_impedance(ml.substrate)

    # Extract transmission and reflection coefficients
    # Same formalism as photonic: t = 1/M[1,1], r = M[2,1]/M[1,1]
    t = 1 / M[1, 1]
    r = M[2, 1] / M[1, 1]

    # Power coefficients
    R = abs2(r)

    # Transmission with impedance correction
    # T = (Z_inc/Z_sub) * |t|² for energy flux
    T = (Z_inc / Z_sub) * abs2(t)

    # Lossless: A = 0 (no absorption for elastic materials)
    A = max(0.0, 1.0 - R - T)

    return TMMResult(R, T, A, r, t)
end
