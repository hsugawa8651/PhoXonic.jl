# Last-Modified: 2025-12-14T18:00:51+09:00

#=
Wave types for photonic and phononic crystals.
=#

"""
Abstract type for all wave types.
"""
abstract type WaveType end

"""
Abstract type for photonic (electromagnetic) waves.
"""
abstract type PhotonicWave <: WaveType end

"""
Abstract type for phononic (elastic) waves.
"""
abstract type PhononicWave <: WaveType end

# ============================================================================
# Photonic Wave Types
# ============================================================================

"""
    TEWave <: PhotonicWave

Transverse Electric polarization (H_z mode) for 2D photonic crystals.

Electric field is in the xy-plane, magnetic field H is along z.
The eigenvalue equation is: ∇×(ε⁻¹∇×H_z) = (ω/c)² μ H_z

# Field Components
- Solved: H_z (scalar)
- Derived: E_x, E_y from H_z

# Example
```julia
solver = Solver(TEWave(), geo, (64, 64))
solver = Solver(TEWave(), geo, (64, 64), KrylovKitMethod())
solver = Solver(TEWave(), geo, (64, 64), LOBPCGMethod())
```

See also: [`TMWave`](@ref), [`FullVectorEM`](@ref)
"""
struct TEWave <: PhotonicWave end

"""
    TMWave <: PhotonicWave

Transverse Magnetic polarization (E_z mode) for 2D photonic crystals.

Magnetic field is in the xy-plane, electric field E is along z.
The eigenvalue equation is: ∇×(μ⁻¹∇×E_z) = (ω/c)² ε E_z

# Field Components
- Solved: E_z (scalar)
- Derived: H_x, H_y from E_z

# Example
```julia
solver = Solver(TMWave(), geo, (64, 64))
solver = Solver(TMWave(), geo, (64, 64), KrylovKitMethod())
solver = Solver(TMWave(), geo, (64, 64), LOBPCGMethod())
```

See also: [`TEWave`](@ref), [`FullVectorEM`](@ref)
"""
struct TMWave <: PhotonicWave end

"""
    Photonic1D <: PhotonicWave

Scalar electromagnetic wave for 1D photonic structures (Bragg reflectors, etc.).

The eigenvalue equation is: -d/dx(ε⁻¹ d/dx E) = (ω/c)² E

# Example
```julia
lat = lattice_1d(1.0)
geo = Geometry(lat, mat1, [(Segment(0.0, 0.5), mat2)])
solver = Solver(Photonic1D(), geo, 128; cutoff=20)
solver = Solver(Photonic1D(), geo, 128, KrylovKitMethod(); cutoff=20)
solver = Solver(Photonic1D(), geo, 128, LOBPCGMethod(); cutoff=20)
```

See also: [`TEWave`](@ref), [`TMWave`](@ref), [`Longitudinal1D`](@ref)
"""
struct Photonic1D <: PhotonicWave end

"""
    FullVectorEM <: PhotonicWave

Full vector electromagnetic wave for 3D photonic crystals.

Uses H-field formulation: ∇×(ε⁻¹∇×H) = (ω/c)² μ H

# Field Components
- Solved: H_x, H_y, H_z (3 components)
- Physical modes: 2 transverse modes per k-point

# Notes
- Produces spurious longitudinal modes at ω ≈ 0 (unphysical)
- Use `shift` parameter with iterative solvers to skip spurious modes

# Example
```julia
lat = cubic_lattice(1.0)
geo = Geometry(lat, air, [(Sphere([0,0,0], 0.3), rod)])
solver = Solver(FullVectorEM(), geo, (16, 16, 16); cutoff=3)

# Iterative solvers require shift to skip spurious modes
solver = Solver(FullVectorEM(), geo, (16, 16, 16), KrylovKitMethod(shift=0.01); cutoff=3)
solver = Solver(FullVectorEM(), geo, (16, 16, 16), LOBPCGMethod(shift=0.01); cutoff=3)
```

See also: [`TEWave`](@ref), [`TMWave`](@ref), [`FullElastic`](@ref)
"""
struct FullVectorEM <: PhotonicWave end

# ============================================================================
# Phononic Wave Types
# ============================================================================

"""
    SHWave <: PhononicWave

Shear Horizontal (out-of-plane, anti-plane) elastic wave for 2D phononic crystals.

Displacement u_z is perpendicular to the xy-plane.
The eigenvalue equation is: ∇·(C₄₄∇u_z) = -ρω² u_z

# Field Components
- Solved: u_z (scalar)
- Decoupled from in-plane (P-SV) modes

# Notes
- Eigenvalues ω² can be O(10¹⁰) for typical materials
- KrylovKitMethod uses automatic scaling for numerical stability
- LOBPCGMethod works without explicit scaling

# Example
```julia
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
geo = Geometry(lat, epoxy, [(Circle([0,0], 0.3), steel)])

solver = Solver(SHWave(), geo, (64, 64))
solver = Solver(SHWave(), geo, (64, 64), KrylovKitMethod())
solver = Solver(SHWave(), geo, (64, 64), LOBPCGMethod())
```

See also: [`PSVWave`](@ref), [`FullElastic`](@ref)
"""
struct SHWave <: PhononicWave end

"""
    PSVWave <: PhononicWave

P-SV (in-plane) elastic wave for 2D phononic crystals.

Displacement (u_x, u_y) is in the xy-plane, coupling P (longitudinal) and
SV (shear vertical) polarizations.

# Field Components
- Solved: u_x, u_y (2 components)
- Produces 2 bands per k-point (quasi-P and quasi-SV)

# Notes
- Eigenvalues ω² can be O(10¹⁰) for typical materials
- KrylovKitMethod uses automatic scaling for numerical stability
- LOBPCGMethod works without explicit scaling

# Example
```julia
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
geo = Geometry(lat, epoxy, [(Circle([0,0], 0.3), steel)])

solver = Solver(PSVWave(), geo, (64, 64))
solver = Solver(PSVWave(), geo, (64, 64), KrylovKitMethod())
solver = Solver(PSVWave(), geo, (64, 64), LOBPCGMethod())
```

See also: [`SHWave`](@ref), [`FullElastic`](@ref)
"""
struct PSVWave <: PhononicWave end

"""
    Longitudinal1D <: PhononicWave

Longitudinal elastic wave for 1D phononic structures (superlattices, etc.).

The eigenvalue equation is: d/dx(C₁₁ du/dx) = -ρω² u

# Notes
- Eigenvalues ω² can be O(10¹⁰) for typical materials
- KrylovKitMethod uses automatic scaling for numerical stability
- LOBPCGMethod works without explicit scaling

# Example
```julia
lat = lattice_1d(1.0)
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
geo = Geometry(lat, epoxy, [(Segment(0.0, 0.5), steel)])

solver = Solver(Longitudinal1D(), geo, 128; cutoff=20)
solver = Solver(Longitudinal1D(), geo, 128, KrylovKitMethod(); cutoff=20)
solver = Solver(Longitudinal1D(), geo, 128, LOBPCGMethod(); cutoff=20)
```

See also: [`SHWave`](@ref), [`PSVWave`](@ref), [`Photonic1D`](@ref)
"""
struct Longitudinal1D <: PhononicWave end

"""
    FullElastic <: PhononicWave

Full elastic wave for 3D phononic crystals.

Solves the elastodynamic equation: ∇·σ = -ρω²u where σ = C:ε

# Field Components
- Solved: u_x, u_y, u_z (3 components)
- Produces 3 bands per k-point (1 quasi-P + 2 quasi-S)

# Notes
- Eigenvalues ω² can be O(10¹⁰) for typical materials
- Use `shift` parameter with iterative solvers for 3D

# Example
```julia
lat = cubic_lattice(1.0)
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
geo = Geometry(lat, epoxy, [(Sphere([0,0,0], 0.3), steel)])

solver = Solver(FullElastic(), geo, (16, 16, 16); cutoff=3)

# Iterative solvers may require shift for 3D
solver = Solver(FullElastic(), geo, (16, 16, 16), KrylovKitMethod(shift=0.01); cutoff=3)
solver = Solver(FullElastic(), geo, (16, 16, 16), LOBPCGMethod(shift=0.01); cutoff=3)
```

See also: [`SHWave`](@ref), [`PSVWave`](@ref), [`FullVectorEM`](@ref)
"""
struct FullElastic <: PhononicWave end

# ============================================================================
# Wave structure traits (for dispatch without Union types)
# ============================================================================

"""
Trait types for wave structure, enabling dispatch based on dimension and components.
"""
struct ScalarWave1D end
struct ScalarWave2D end
struct VectorWave2D end   # 2 components (PSVWave)
struct VectorWave3D end   # 3 components (FullVectorEM, FullElastic)

"""
    wave_structure(::WaveType)

Return the structure trait for a wave type.
Enables dispatch based on dimensionality and number of components.
"""
# Error fallback for unsupported wave types
function wave_structure(w::WaveType)
    throw(ArgumentError("Unsupported wave type for wave_structure: $(typeof(w))"))
end

wave_structure(::TEWave) = ScalarWave2D()
wave_structure(::TMWave) = ScalarWave2D()
wave_structure(::SHWave) = ScalarWave2D()
wave_structure(::PSVWave) = VectorWave2D()
wave_structure(::FullVectorEM) = VectorWave3D()
wave_structure(::FullElastic) = VectorWave3D()
wave_structure(::Photonic1D) = ScalarWave1D()
wave_structure(::Longitudinal1D) = ScalarWave1D()

# ============================================================================
# Wave properties
# ============================================================================

"""
    ncomponents(::WaveType)

Return the number of field components for a wave type.
"""
# Error fallback for unsupported wave types
function ncomponents(w::WaveType)
    throw(ArgumentError("Unsupported wave type for ncomponents: $(typeof(w))"))
end

ncomponents(::TEWave) = 1
ncomponents(::TMWave) = 1
ncomponents(::Photonic1D) = 1
ncomponents(::SHWave) = 1
ncomponents(::Longitudinal1D) = 1
ncomponents(::PSVWave) = 2
ncomponents(::FullVectorEM) = 3
ncomponents(::FullElastic) = 3

"""
    applicable_dimension(::Type{W}) where W <: WaveType

Return the applicable dimension type for a wave type.
"""
# Error fallback for unsupported wave types
function applicable_dimension(::Type{W}) where {W<:WaveType}
    throw(ArgumentError("Unsupported wave type for applicable_dimension: $W"))
end

applicable_dimension(::Type{TEWave}) = Dim2
applicable_dimension(::Type{TMWave}) = Dim2
applicable_dimension(::Type{SHWave}) = Dim2
applicable_dimension(::Type{PSVWave}) = Dim2
applicable_dimension(::Type{Photonic1D}) = Dim1
applicable_dimension(::Type{Longitudinal1D}) = Dim1
applicable_dimension(::Type{FullVectorEM}) = Dim3
applicable_dimension(::Type{FullElastic}) = Dim3

"""
    is_scalar(::WaveType)

Return true if the wave type has a scalar (single component) field.
"""
is_scalar(w::WaveType) = ncomponents(w) == 1
