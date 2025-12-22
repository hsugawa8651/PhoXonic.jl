# Transfer matrix calculations for normal and oblique incidence

"""
    fresnel_coefficients(n1, n2) -> (r, t)

Calculate Fresnel reflection and transmission coefficients for
normal incidence at an interface between media with refractive
indices n1 (incident) and n2 (transmitted).

For electric field amplitudes:
- r = (n1 - n2) / (n1 + n2)
- t = 2n1 / (n1 + n2)

Power reflectance and transmittance:
- R = |r|²
- T = (n2/n1) |t|²
- R + T = 1

Both real and complex refractive indices are supported.
"""
function fresnel_coefficients(n1::Number, n2::Number)
    r = (n1 - n2) / (n1 + n2)
    t = 2 * n1 / (n1 + n2)
    return r, t
end

"""
    propagation_matrix(n, d, λ) -> Matrix{ComplexF64}

Calculate the 2×2 propagation matrix for a layer with refractive index n
and thickness d at wavelength λ.

The propagation matrix accounts for phase accumulation:
    P = [exp(iδ)  0    ]
        [0        exp(-iδ)]

where δ = 2πnd/λ is the optical phase.

For negative d (backward propagation), exp(-iδ) becomes exp(iδ),
so P(-d) = P(d)⁻¹.

Complex n is supported for absorbing materials.
"""
function propagation_matrix(n::Number, d::Real, λ::Real)
    δ = 2π * n * d / λ
    return ComplexF64[exp(im * δ) 0; 0 exp(-im * δ)]
end

"""
    interface_matrix(n1, n2) -> Matrix{ComplexF64}

Calculate the 2×2 interface matrix for light going from medium n1
to medium n2 at normal incidence.

The interface matrix is:
    M = (1/t) * [1  r]
                [r  1]

where r and t are the Fresnel reflection and transmission coefficients.

Properties:
- det(M) = n2/n1
- M is symmetric

Both real and complex refractive indices are supported.
"""
function interface_matrix(n1::Number, n2::Number)
    r, t = fresnel_coefficients(n1, n2)
    return ComplexF64[1 r; r 1] / t
end

"""
    system_matrix(ml::Multilayer, λ::Real) -> Matrix{ComplexF64}

Calculate the total 2×2 transfer matrix for a multilayer structure
at wavelength λ for normal incidence.

The system matrix M relates field amplitudes:
    [E_inc^+]     [E_sub^+]
    [E_inc^-] = M [E_sub^-]

where E^+ and E^- are forward and backward propagating waves.

The matrix is computed as a product from incident to substrate:
    M = I_{0,1} × P_1 × I_{1,2} × P_2 × ... × I_{N,N+1}

where I_{i,j} are interface matrices and P_j are propagation matrices.
"""
function system_matrix(ml::Multilayer, λ::Real)
    # Get refractive indices
    n_inc = _get_n(ml.incident)
    n_sub = _get_n(ml.substrate)

    # Start with identity
    M = ComplexF64[1 0; 0 1]

    # Previous refractive index (starting from incident medium)
    n_prev = n_inc

    # Loop through layers
    for layer in ml.layers
        n_layer = _get_n(layer.material)
        d = thickness(layer)

        # Interface from previous medium to this layer
        M = M * interface_matrix(n_prev, n_layer)

        # Propagation through this layer
        M = M * propagation_matrix(n_layer, d, λ)

        n_prev = n_layer
    end

    # Final interface to substrate
    M = M * interface_matrix(n_prev, n_sub)

    return M
end

# Helper function to get refractive index from material
function _get_n(mat::Dielectric)
    # n = √ε for non-magnetic materials
    sqrt(mat.ε)
end

function _get_n(mat::LossyDielectric)
    # n = √(εμ) for complex materials
    sqrt(mat.ε * mat.μ)
end

# ============================================================================
# Phononic (elastic wave) support
# ============================================================================

"""
    acoustic_impedance(mat::IsotropicElastic) -> Float64

Calculate the acoustic impedance Z = ρ * c_L for longitudinal waves.

# Arguments
- `mat`: Isotropic elastic material

# Returns
- `Z`: Acoustic impedance [kg/(m²·s)] = [Rayl]
"""
function acoustic_impedance(mat::IsotropicElastic)
    ρ = density(mat)
    c = longitudinal_velocity(mat)
    return ρ * c
end

function acoustic_impedance(mat::ElasticVoid)
    ρ = density(mat)
    c = longitudinal_velocity(mat)
    return ρ * c
end

"""
    acoustic_fresnel(Z1, Z2) -> (r, t)

Calculate acoustic reflection and transmission coefficients at an interface.

For stress/pressure amplitude:
- r = (Z2 - Z1) / (Z2 + Z1)
- t = 2*Z2 / (Z1 + Z2)

For intensity (power):
- R = |r|²
- T = 4*Z1*Z2 / (Z1 + Z2)²
- R + T = 1

# Arguments
- `Z1`: Acoustic impedance of incident medium
- `Z2`: Acoustic impedance of transmitted medium

# Returns
- `r`: Amplitude reflection coefficient
- `t`: Amplitude transmission coefficient
"""
function acoustic_fresnel(Z1::Real, Z2::Real)
    r = (Z2 - Z1) / (Z2 + Z1)
    t = 2 * Z2 / (Z1 + Z2)
    return r, t
end

"""
    acoustic_propagation_matrix(mat, d, λ) -> Matrix{ComplexF64}

Calculate the 2×2 propagation matrix for an elastic layer.

# Arguments
- `mat`: Elastic material (IsotropicElastic or ElasticVoid)
- `d`: Layer thickness
- `λ`: Wavelength (in physical units matching d)
"""
function acoustic_propagation_matrix(
    mat::Union{IsotropicElastic,ElasticVoid}, d::Real, λ::Real
)
    c = longitudinal_velocity(mat)
    # Phase: δ = k*d = 2π*d/λ (wavelength in the material)
    # For physical wavelength λ in incident medium, we need λ_mat = λ * c_mat / c_inc
    # But for TMM we use wavelength directly
    δ = 2π * d / λ
    return ComplexF64[exp(im * δ) 0; 0 exp(-im * δ)]
end

"""
    acoustic_interface_matrix(mat1, mat2) -> Matrix{ComplexF64}

Calculate the 2×2 interface matrix for elastic waves.

# Arguments
- `mat1`: Material on incident side
- `mat2`: Material on transmitted side
"""
function acoustic_interface_matrix(
    mat1::Union{IsotropicElastic,ElasticVoid}, mat2::Union{IsotropicElastic,ElasticVoid}
)
    Z1 = acoustic_impedance(mat1)
    Z2 = acoustic_impedance(mat2)
    r, t = acoustic_fresnel(Z1, Z2)
    return ComplexF64[1 r; r 1] / t
end

"""
    system_matrix_acoustic(ml::Multilayer, λ::Real)

Calculate the total transfer matrix for a phononic multilayer structure.

# Arguments
- `ml`: Multilayer structure with elastic materials
- `λ`: Wavelength (physical units)

# Returns
- 2×2 complex transfer matrix
"""
function system_matrix_acoustic(ml::Multilayer{<:ElasticMaterial}, λ::Real)
    # Start with identity
    M = ComplexF64[1 0; 0 1]

    # Previous material
    mat_prev = ml.incident

    # Loop through layers
    for layer in ml.layers
        mat_layer = layer.material
        d = thickness(layer)

        # Interface from previous to layer
        M = M * acoustic_interface_matrix(mat_prev, mat_layer)

        # Propagation through layer
        M = M * acoustic_propagation_matrix(mat_layer, d, λ)

        mat_prev = mat_layer
    end

    # Final interface to substrate
    M = M * acoustic_interface_matrix(mat_prev, ml.substrate)

    return M
end

# ============================================================================
# Oblique incidence
# ============================================================================

"""
    snell_angle(n1, n2, θ1) -> θ2

Calculate the refraction angle θ2 using Snell's law: n1 sin(θ1) = n2 sin(θ2).

Returns a complex angle for total internal reflection (when sin(θ2) > 1).

# Arguments
- `n1`: Refractive index of incident medium
- `n2`: Refractive index of transmitted medium
- `θ1`: Incident angle in radians

# Returns
- `θ2`: Transmitted angle in radians (may be complex for TIR)
"""
function snell_angle(n1::Number, n2::Number, θ1::Number)
    sin_θ2 = n1 * sin(θ1) / n2
    # Convert to complex to handle total internal reflection (|sin_θ2| > 1)
    return asin(ComplexF64(sin_θ2))
end

"""
    fresnel_oblique(n1, n2, θ1, polarization) -> (r, t)

Calculate Fresnel reflection and transmission coefficients for
oblique incidence at an interface.

# Arguments
- `n1`: Refractive index of incident medium
- `n2`: Refractive index of transmitted medium
- `θ1`: Incident angle in radians
- `polarization`: `:TE` (s-wave, E perpendicular to plane) or `:TM` (p-wave, E in plane)

# Returns
- `r`: Complex reflection coefficient
- `t`: Complex transmission coefficient

# Formulas
TE (s-wave):
- r_s = (n1 cos θ1 - n2 cos θ2) / (n1 cos θ1 + n2 cos θ2)
- t_s = 2 n1 cos θ1 / (n1 cos θ1 + n2 cos θ2)

TM (p-wave):
- r_p = (n2 cos θ1 - n1 cos θ2) / (n2 cos θ1 + n1 cos θ2)
- t_p = 2 n1 cos θ1 / (n2 cos θ1 + n1 cos θ2)
"""
function fresnel_oblique(n1::Number, n2::Number, θ1::Number, polarization::Symbol)
    θ2 = snell_angle(n1, n2, θ1)
    cos_θ1 = cos(θ1)
    cos_θ2 = cos(θ2)

    if polarization == :TE
        # TE (s-wave): E-field perpendicular to plane of incidence
        r = (n1 * cos_θ1 - n2 * cos_θ2) / (n1 * cos_θ1 + n2 * cos_θ2)
        t = 2 * n1 * cos_θ1 / (n1 * cos_θ1 + n2 * cos_θ2)
    elseif polarization == :TM
        # TM (p-wave): E-field in plane of incidence
        r = (n2 * cos_θ1 - n1 * cos_θ2) / (n2 * cos_θ1 + n1 * cos_θ2)
        t = 2 * n1 * cos_θ1 / (n2 * cos_θ1 + n1 * cos_θ2)
    else
        error("Invalid polarization: $polarization. Use :TE or :TM")
    end

    return r, t
end

"""
    propagation_matrix_oblique(n, d, λ, θ) -> Matrix{ComplexF64}

Calculate the 2×2 propagation matrix for a layer at oblique incidence.

The phase is modified by the component of the wave vector perpendicular
to the interface: δ = 2π n d cos(θ) / λ

# Arguments
- `n`: Refractive index of the layer
- `d`: Layer thickness
- `λ`: Wavelength
- `θ`: Angle in the layer (from Snell's law)
"""
function propagation_matrix_oblique(n::Number, d::Real, λ::Real, θ::Number)
    # Phase accumulation perpendicular to interface
    δ = 2π * n * d * cos(θ) / λ
    return ComplexF64[exp(im * δ) 0; 0 exp(-im * δ)]
end

"""
    interface_matrix_oblique(n1, n2, θ1, polarization) -> Matrix{ComplexF64}

Calculate the 2×2 interface matrix for oblique incidence.

# Arguments
- `n1`: Refractive index of incident medium
- `n2`: Refractive index of transmitted medium
- `θ1`: Incident angle
- `polarization`: `:TE` or `:TM`
"""
function interface_matrix_oblique(n1::Number, n2::Number, θ1::Number, polarization::Symbol)
    r, t = fresnel_oblique(n1, n2, θ1, polarization)
    return ComplexF64[1 r; r 1] / t
end

"""
    system_matrix_oblique(ml::Multilayer, λ::Real, θ_inc::Real, polarization::Symbol)

Calculate the total transfer matrix for a multilayer at oblique incidence.

# Arguments
- `ml`: Multilayer structure
- `λ`: Wavelength
- `θ_inc`: Incident angle in the first medium (radians)
- `polarization`: `:TE` or `:TM`

# Returns
- 2×2 complex transfer matrix
"""
function system_matrix_oblique(ml::Multilayer, λ::Real, θ_inc::Real, polarization::Symbol)
    # Get refractive indices
    n_inc = _get_n(ml.incident)
    n_sub = _get_n(ml.substrate)

    # Start with identity
    M = ComplexF64[1 0; 0 1]

    # Track angle and refractive index through layers
    n_prev = n_inc
    θ_prev = θ_inc

    # Loop through layers
    for layer in ml.layers
        n_layer = _get_n(layer.material)
        d = thickness(layer)

        # Angle in this layer (Snell's law)
        θ_layer = snell_angle(n_prev, n_layer, θ_prev)

        # Interface from previous medium to this layer
        M = M * interface_matrix_oblique(n_prev, n_layer, θ_prev, polarization)

        # Propagation through this layer
        M = M * propagation_matrix_oblique(n_layer, d, λ, θ_layer)

        n_prev = n_layer
        θ_prev = θ_layer
    end

    # Final interface to substrate
    M = M * interface_matrix_oblique(n_prev, n_sub, θ_prev, polarization)

    return M
end
