# Last-Modified: 2025-12-09T22:36:35+09:00

#=
Convolution matrix construction from real-space data.
=#

"""
    convolution_matrix(f::AbstractMatrix, basis::PlaneWaveBasis{Dim2})

Compute the convolution matrix for a 2D function f(r).

The convolution matrix C has elements:
    C[i,j] = f̂(G_i - G_j)
where f̂ is the Fourier transform of f.

# Arguments
- `f`: Real-space function values on a grid (Nx × Ny)
- `basis`: Plane wave basis

# Returns
A complex matrix of size (num_pw × num_pw).
"""
function convolution_matrix(f::AbstractMatrix, basis::PlaneWaveBasis{Dim2})
    Nx, Ny = size(f)
    N = basis.num_pw

    # Compute FFT of f (with proper normalization and shifting)
    # fftshift moves zero-frequency to center
    f_fft = fftshift(fft(fftshift(f))) / length(f)

    # Center indices for accessing FFT result
    cx = (Nx + 1) ÷ 2 + 1
    cy = (Ny + 1) ÷ 2 + 1

    # Build convolution matrix
    C = zeros(ComplexF64, N, N)

    for i in 1:N
        pi, qi = basis.indices[i]
        for j in 1:N
            pj, qj = basis.indices[j]
            # Difference in indices
            dp = pi - pj
            dq = qi - qj
            # Check bounds
            idx_x = cx + dp
            idx_y = cy + dq
            if 1 <= idx_x <= Nx && 1 <= idx_y <= Ny
                C[i, j] = f_fft[idx_x, idx_y]
            end
        end
    end

    return C
end

"""
    convolution_matrix(f::AbstractVector, basis::PlaneWaveBasis{Dim1})

Compute the convolution matrix for a 1D function f(x).
"""
function convolution_matrix(f::AbstractVector, basis::PlaneWaveBasis{Dim1})
    Nx = length(f)
    N = basis.num_pw

    # Compute FFT of f
    f_fft = fftshift(fft(fftshift(f))) / Nx

    # Center index
    cx = (Nx + 1) ÷ 2 + 1

    # Build convolution matrix
    C = zeros(ComplexF64, N, N)

    for i in 1:N
        pi = basis.indices[i][1]
        for j in 1:N
            pj = basis.indices[j][1]
            dp = pi - pj
            idx = cx + dp
            if 1 <= idx <= Nx
                C[i, j] = f_fft[idx]
            end
        end
    end

    return C
end

"""
    convolution_matrix(f::AbstractArray{T,3}, basis::PlaneWaveBasis{Dim3}) where T

Compute the convolution matrix for a 3D function f(r).
"""
function convolution_matrix(f::AbstractArray{T,3}, basis::PlaneWaveBasis{Dim3}) where {T}
    Nx, Ny, Nz = size(f)
    N = basis.num_pw

    # Compute FFT of f
    f_fft = fftshift(fft(fftshift(f))) / length(f)

    # Center indices
    cx = (Nx + 1) ÷ 2 + 1
    cy = (Ny + 1) ÷ 2 + 1
    cz = (Nz + 1) ÷ 2 + 1

    # Build convolution matrix
    C = zeros(ComplexF64, N, N)

    for i in 1:N
        pi, qi, ri = basis.indices[i]
        for j in 1:N
            pj, qj, rj = basis.indices[j]
            dp = pi - pj
            dq = qi - qj
            dr = ri - rj
            idx_x = cx + dp
            idx_y = cy + dq
            idx_z = cz + dr
            if 1 <= idx_x <= Nx && 1 <= idx_y <= Ny && 1 <= idx_z <= Nz
                C[i, j] = f_fft[idx_x, idx_y, idx_z]
            end
        end
    end

    return C
end

"""
    inverse_convolution_matrix(f::AbstractArray, basis::PlaneWaveBasis)

Compute the convolution matrix of 1/f.
Useful for TM modes where we need ε⁻¹.
"""
function inverse_convolution_matrix(f::AbstractArray, basis::PlaneWaveBasis)
    convolution_matrix(1.0 ./ f, basis)
end
