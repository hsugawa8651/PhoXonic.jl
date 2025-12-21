# Transfer matrix calculations for normal incidence

"""
    fresnel_coefficients(n1::Real, n2::Real) -> (r, t)

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
"""
function fresnel_coefficients(n1::Real, n2::Real)
    r = (n1 - n2) / (n1 + n2)
    t = 2 * n1 / (n1 + n2)
    return r, t
end

"""
    propagation_matrix(n::Real, d::Real, λ::Real) -> Matrix{ComplexF64}

Calculate the 2×2 propagation matrix for a layer with refractive index n
and thickness d at wavelength λ.

The propagation matrix accounts for phase accumulation:
    P = [exp(iδ)  0    ]
        [0        exp(-iδ)]

where δ = 2πnd/λ is the optical phase.

For negative d (backward propagation), exp(-iδ) becomes exp(iδ),
so P(-d) = P(d)⁻¹.
"""
function propagation_matrix(n::Real, d::Real, λ::Real)
    δ = 2π * n * d / λ
    return ComplexF64[exp(im * δ) 0; 0 exp(-im * δ)]
end

"""
    interface_matrix(n1::Real, n2::Real) -> Matrix{ComplexF64}

Calculate the 2×2 interface matrix for light going from medium n1
to medium n2 at normal incidence.

The interface matrix is:
    M = (1/t) * [1  r]
                [r  1]

where r and t are the Fresnel reflection and transmission coefficients.

Properties:
- det(M) = n2/n1
- M is symmetric
"""
function interface_matrix(n1::Real, n2::Real)
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
