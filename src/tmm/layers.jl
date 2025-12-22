# Layer and Multilayer types for Transfer Matrix Method

"""
    Layer{M<:Material}

A single layer in a multilayer structure.

# Fields
- `material::M`: Material of the layer
- `thickness::Float64`: Physical thickness of the layer

# Examples
```julia
glass = Dielectric(2.25)
layer = Layer(glass, 100.0)  # 100 nm thick glass layer
```
"""
struct Layer{M<:Material}
    material::M
    thickness::Float64
end

Layer(mat::Material, d::Real) = Layer(mat, Float64(d))

"""
    thickness(layer::Layer)

Return the physical thickness of the layer.
"""
thickness(layer::Layer) = layer.thickness

"""
    material(layer::Layer)

Return the material of the layer.
"""
material(layer::Layer) = layer.material

"""
    Multilayer{M<:Material}

A multilayer structure consisting of multiple layers between
semi-infinite incident and substrate media.

# Fields
- `layers::Vector{Layer{M}}`: Array of layers
- `incident::M`: Incident medium (semi-infinite)
- `substrate::M`: Substrate medium (semi-infinite)

# Examples
```julia
air = Dielectric(1.0)
glass = Dielectric(2.25)
layer = Layer(glass, 100.0)
ml = Multilayer([layer], air, air)
```
"""
struct Multilayer{M<:Material}
    layers::Vector{Layer{M}}
    incident::M
    substrate::M
end

# Constructor with automatic type promotion
function Multilayer(layers::Vector{<:Layer}, incident::Material, substrate::Material)
    # Promote all materials to common type
    M = promote_type(typeof(incident), typeof(substrate))
    for l in layers
        M = promote_type(M, typeof(material(l)))
    end

    # Convert layers
    converted_layers = [Layer(convert(M, material(l)), thickness(l)) for l in layers]

    Multilayer{M}(converted_layers, convert(M, incident), convert(M, substrate))
end

"""
    layers(ml::Multilayer)

Return the array of layers.
"""
layers(ml::Multilayer) = ml.layers

"""
    incident(ml::Multilayer)

Return the incident medium.
"""
incident(ml::Multilayer) = ml.incident

"""
    substrate(ml::Multilayer)

Return the substrate medium.
"""
substrate(ml::Multilayer) = ml.substrate

"""
    nlayers(ml::Multilayer)

Return the number of layers.
"""
nlayers(ml::Multilayer) = length(ml.layers)

"""
    total_thickness(ml::Multilayer)

Return the total thickness of all layers.
"""
total_thickness(ml::Multilayer) = sum(thickness(l) for l in ml.layers)

# ============================================================================
# Helper functions
# ============================================================================

"""
    periodic_multilayer(unit_cell::Vector{<:Layer}, N::Int;
                        incident=nothing, substrate=nothing)

Create a multilayer by repeating a unit cell N times.

# Arguments
- `unit_cell`: Vector of layers forming one period
- `N`: Number of periods

# Keyword Arguments
- `incident`: Incident medium (defaults to first layer's material type with ε=1 or air-like)
- `substrate`: Substrate medium (defaults to same as incident)
"""
function periodic_multilayer(
    unit_cell::Vector{<:Layer}, N::Int; incident=nothing, substrate=nothing
)
    repeated_layers = repeat(unit_cell, N)

    # Default incident/substrate
    if incident === nothing
        mat = material(unit_cell[1])
        if mat isa Dielectric
            incident = Dielectric(1.0)
        elseif mat isa IsotropicElastic
            # Use first material as default for elastic
            incident = mat
        else
            incident = mat
        end
    end

    if substrate === nothing
        substrate = incident
    end

    Multilayer(repeated_layers, incident, substrate)
end

"""
    bragg_mirror(n_hi::Real, n_lo::Real, λ0::Real, N::Int;
                 incident=Dielectric(1.0), substrate=Dielectric(1.0))

Create a Bragg mirror (quarter-wave stack) for center wavelength λ0.

# Arguments
- `n_hi`: High refractive index
- `n_lo`: Low refractive index
- `λ0`: Center wavelength (same units as returned thicknesses)
- `N`: Number of layer pairs

# Returns
A Multilayer with 2N layers alternating between high and low index materials.
Layer thicknesses satisfy the quarter-wave condition: n*d = λ0/4
"""
function bragg_mirror(
    n_hi::Real,
    n_lo::Real,
    λ0::Real,
    N::Int;
    incident=Dielectric(1.0),
    substrate=Dielectric(1.0),
)
    # Quarter-wave thicknesses
    d_hi = λ0 / (4 * n_hi)
    d_lo = λ0 / (4 * n_lo)

    # Materials (using ε = n²)
    mat_hi = Dielectric(n_hi^2)
    mat_lo = Dielectric(n_lo^2)

    # Unit cell: high-index then low-index
    unit_cell = [Layer(mat_hi, d_hi), Layer(mat_lo, d_lo)]

    periodic_multilayer(unit_cell, N; incident=incident, substrate=substrate)
end
