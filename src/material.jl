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
# LossyDielectric (TMM only)
# ============================================================================

"""
    LossyDielectric(ε, μ=1.0)

Lossy dielectric material with complex permittivity. TMM only.

Complex permittivity ε = ε' + iε'' supports absorption.
Relation to refractive index: ε = n² where n = n' + in''.

# Examples
```julia
# Metallic absorption
gold = LossyDielectric(-10.0 + 1.0im)

# Weak absorption
glass_absorbing = LossyDielectric(2.25 + 0.01im)
```

# Note
Cannot be used with PWE (band structure). PWE requires Hermitian eigenvalue
problems, which need real ε. Complex ε produces complex eigenvalues.
"""
struct LossyDielectric <: PhotonicMaterial
    ε::ComplexF64
    μ::ComplexF64
end

LossyDielectric(ε::Number) = LossyDielectric(ComplexF64(ε), ComplexF64(1.0))
LossyDielectric(ε::Number, μ::Number) = LossyDielectric(ComplexF64(ε), ComplexF64(μ))

"""
    permittivity(m::LossyDielectric)

Return the complex permittivity.
"""
permittivity(m::LossyDielectric) = m.ε

"""
    permeability(m::LossyDielectric)

Return the complex permeability.
"""
permeability(m::LossyDielectric) = m.μ

"""
    refractive_index(m::LossyDielectric)

Return the complex refractive index n = √(εμ).
"""
refractive_index(m::LossyDielectric) = sqrt(m.ε * m.μ)

# Type conversion: Dielectric → LossyDielectric
function Base.convert(::Type{LossyDielectric}, d::Dielectric)
    LossyDielectric(ComplexF64(d.ε), ComplexF64(d.μ))
end

# Type promotion rules
Base.promote_rule(::Type{Dielectric}, ::Type{LossyDielectric}) = LossyDielectric

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
    IsotropicElastic(; ρ=ρ, λ=λ, μ=μ)
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

# ============================================================================
# ElasticVoid (Tanaka Limit)
# ============================================================================

"""
    ElasticVoid(; ρ_ratio=1e-7)

Void region for phononic crystals using the Tanaka limit approach.

The Tanaka limit (Tanaka et al., Phys. Rev. B 62, 7387 (2000)) handles void
regions in PWE by setting ρ/C → 0, which pushes spurious flat bands to
high frequency (ω ~ √(C/ρ) → ∞).

The implementation uses normalized elastic constants (C = 1) and small
density (ρ = ρ_ratio), giving transverse velocity v_T = √(μ/ρ) = √(1/ρ_ratio).

# Arguments
- `ρ_ratio::Real=1e-7`: Density ratio. Smaller values push spurious bands
  higher but may affect numerical stability. Default 1e-7 gives v_T ≈ 3162.

# Examples
```julia
void = ElasticVoid()              # Default: ρ = 1e-7, v_T ≈ 3162
void = ElasticVoid(ρ_ratio=1e-8)  # Higher v_T = 10000
```

# References
- Tanaka et al., Phys. Rev. B 62, 7387 (2000)
- Maldovan & Thomas, Appl. Phys. B 83, 595 (2006)
"""
struct ElasticVoid <: ElasticMaterial
    ρ::Float64    # Density (very small)
    C11::Float64  # λ + 2μ
    C12::Float64  # λ
    C44::Float64  # μ (shear modulus)
end

function ElasticVoid(; ρ_ratio::Real=1e-7)
    # Normalized elastic constants: μ = 1, λ = 1
    μ = 1.0
    λ = 1.0
    ρ = ρ_ratio * μ  # Small density for Tanaka limit
    ElasticVoid(Float64(ρ), Float64(λ + 2μ), Float64(λ), Float64(μ))
end

"""
    density(m::ElasticVoid)

Return the mass density of the void material.
"""
density(m::ElasticVoid) = m.ρ

"""
    shear_modulus(m::ElasticVoid)

Return the shear modulus (μ = C44).
"""
shear_modulus(m::ElasticVoid) = m.C44

"""
    longitudinal_velocity(m::ElasticVoid)

Return the longitudinal wave velocity v_L = √(C11/ρ).
"""
longitudinal_velocity(m::ElasticVoid) = sqrt(m.C11 / m.ρ)

"""
    transverse_velocity(m::ElasticVoid)

Return the transverse (shear) wave velocity v_T = √(C44/ρ).
For Tanaka limit: v_T → ∞ as ρ → 0.
"""
transverse_velocity(m::ElasticVoid) = sqrt(m.C44 / m.ρ)

# ============================================================================
# Type Promotion for Mixed Elastic Materials
# ============================================================================

# Enable Geometry to mix IsotropicElastic and ElasticVoid
Base.promote_rule(::Type{IsotropicElastic}, ::Type{ElasticVoid}) = ElasticMaterial
Base.promote_rule(::Type{ElasticVoid}, ::Type{IsotropicElastic}) = ElasticMaterial

# Convert to common supertype
Base.convert(::Type{ElasticMaterial}, x::IsotropicElastic) = x
Base.convert(::Type{ElasticMaterial}, x::ElasticVoid) = x
