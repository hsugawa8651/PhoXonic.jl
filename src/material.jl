# Last-Modified: 2025-12-09T22:33:02+09:00

#=
Material definitions for photonic and phononic crystals.
=#

"""
Abstract type for all materials.
"""
abstract type Material end

"""
Abstract type for photonic (electromagnetic) materials.
"""
abstract type PhotonicMaterial <: Material end

"""
Abstract type for elastic (phononic) materials.
"""
abstract type ElasticMaterial <: Material end

# ============================================================================
# Photonic Materials
# ============================================================================

"""
    Dielectric(ε, μ=1.0)

Isotropic dielectric material with permittivity `ε` and permeability `μ`.

# Examples
```julia
air = Dielectric(1.0)
silicon = Dielectric(11.7)
```
"""
struct Dielectric <: PhotonicMaterial
    ε::Float64  # Relative permittivity
    μ::Float64  # Relative permeability
end

Dielectric(ε::Real) = Dielectric(Float64(ε), 1.0)

"""
    permittivity(m::Dielectric)

Return the relative permittivity of the material.
"""
permittivity(m::Dielectric) = m.ε

"""
    permeability(m::Dielectric)

Return the relative permeability of the material.
"""
permeability(m::Dielectric) = m.μ

"""
    refractive_index(m::Dielectric)

Return the refractive index n = √(εμ).
"""
refractive_index(m::Dielectric) = sqrt(m.ε * m.μ)

# ============================================================================
# Elastic Materials
# ============================================================================

"""
    IsotropicElastic(ρ, C11, C12, C44)

Isotropic elastic material with density `ρ` and elastic constants.

For isotropic materials: C44 = (C11 - C12) / 2

# Fields
- `ρ`: Mass density [kg/m³]
- `C11`: Elastic constant (λ + 2μ) [Pa]
- `C12`: Elastic constant (λ) [Pa]
- `C44`: Shear modulus (μ) [Pa]

# Examples
```julia
# Steel
steel = IsotropicElastic(7800.0, 282e9, 113e9, 84.5e9)

# From Lamé parameters
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
```
"""
struct IsotropicElastic <: ElasticMaterial
    ρ::Float64    # Density
    C11::Float64  # λ + 2μ
    C12::Float64  # λ
    C44::Float64  # μ (shear modulus)
end

"""
    IsotropicElastic(; ρ, λ, μ)

Create isotropic elastic material from Lamé parameters.
"""
function IsotropicElastic(; ρ::Real, λ::Real, μ::Real)
    IsotropicElastic(Float64(ρ), Float64(λ + 2μ), Float64(λ), Float64(μ))
end

"""
    from_E_ν(ρ, E, ν)

Create isotropic elastic material from Young's modulus `E` and Poisson's ratio `ν`.
"""
function from_E_ν(ρ::Real, E::Real, ν::Real)
    λ = E * ν / ((1 + ν) * (1 - 2ν))
    μ = E / (2 * (1 + ν))
    IsotropicElastic(ρ=ρ, λ=λ, μ=μ)
end

"""
    density(m::IsotropicElastic)

Return the mass density.
"""
density(m::IsotropicElastic) = m.ρ

"""
    shear_modulus(m::IsotropicElastic)

Return the shear modulus (μ = C44).
"""
shear_modulus(m::IsotropicElastic) = m.C44

"""
    longitudinal_velocity(m::IsotropicElastic)

Return the longitudinal wave velocity v_L = √(C11/ρ).
"""
longitudinal_velocity(m::IsotropicElastic) = sqrt(m.C11 / m.ρ)

"""
    transverse_velocity(m::IsotropicElastic)

Return the transverse (shear) wave velocity v_T = √(C44/ρ).
"""
transverse_velocity(m::IsotropicElastic) = sqrt(m.C44 / m.ρ)
