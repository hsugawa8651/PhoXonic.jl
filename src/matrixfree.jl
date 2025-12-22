# Last-Modified: 2025-12-13T17:25:23+09:00

#=
Matrix-free operators for iterative solvers.

The Hamiltonian H is applied without explicitly constructing the full matrix.
This reduces memory from O(N²) to O(N) and computation from O(N²) to O(N log N).

For photonic/phononic crystals, the eigenvalue problem is:
    LHS * ψ = ω² * RHS * ψ

where:
    LHS = K * M * K  (with K = k + G diagonal, M = convolution matrix)
    RHS = convolution matrix of ρ or μ

The matrix-free approach computes H*v = RHS⁻¹ * LHS * v via FFT:
    1. Apply K in Fourier space (diagonal, O(N))
    2. FFT to real space (O(N log N))
    3. Multiply by material property in real space (O(N))
    4. FFT back to Fourier space (O(N log N))
    5. Apply K in Fourier space (diagonal, O(N))
    6. Solve RHS * y = result (or apply RHS⁻¹)

## Thread Safety (FFTW)

FFTW has specific thread-safety constraints:
- **Plan creation (plan_fft!, etc.)**: NOT thread-safe
- **Plan execution (fft!, mul!, etc.)**: Thread-safe

The design separates concerns for safe parallel usage:
- `FFTContext`: Holds FFT plans. Create ONCE from a single thread before
  entering parallel regions. Can be shared across threads for execution.
- `MatrixFreeWorkspace`: Holds work arrays (input/output buffers).
  Each thread MUST have its own workspace to avoid data races.
- `MatrixFreeOperator`: Combines solver, k-vector, context, and workspace.

Example for parallel k-point loop:
```julia
ctx = FFTContext(solver)  # Single thread, before @threads
Threads.@threads for k in k_points
    ws = MatrixFreeWorkspace(ctx)  # Each thread creates its own
    op = MatrixFreeOperator(solver, k, ctx, ws)
    # ... use op ...
end
```

References:
- https://www.fftw.org/fftw3_doc/Thread-safety.html
- https://github.com/JuliaLang/julia/issues/17972
=#

using LinearMaps

# ============================================================================
# FFT Context and Workspace (for thread-safe FFT plan reuse)
# ============================================================================

"""
    FFTContext{N, T<:Complex, F, I}

Thread-safe container for FFT plans. Plans can be shared across threads
since `fftw_execute` is thread-safe. Create once and reuse.

# Fields
- `resolution`: Grid resolution tuple
- `fft_plan`: Pre-computed forward FFT plan
- `ifft_plan`: Pre-computed inverse FFT plan
"""
struct FFTContext{N,T<:Complex,F,I}
    resolution::NTuple{N,Int}
    fft_plan::F
    ifft_plan::I
end

"""
    FFTContext(resolution::NTuple{N, Int}, ::Type{T}=ComplexF64) where {N, T}

Create an FFT context for the given resolution. Call from a single thread.
"""
function FFTContext(resolution::NTuple{N,Int}, (::Type{T})=ComplexF64) where {N,T<:Complex}
    work = zeros(T, resolution)
    fft_plan = plan_fft!(work)
    ifft_plan = plan_ifft!(work)
    FFTContext{N,T,typeof(fft_plan),typeof(ifft_plan)}(resolution, fft_plan, ifft_plan)
end

"""
    FFTContext(solver::Solver)

Create an FFT context from a solver's resolution.
"""
FFTContext(solver::Solver) = FFTContext(solver.resolution, ComplexF64)

"""
    MatrixFreeWorkspace{N, T<:Complex}

Thread-local workspace arrays for matrix-free operations.
Each thread should have its own workspace instance.

# Fields
- `work_real`: Workspace array in real space
- `work_fourier`: Workspace array in Fourier space
"""
struct MatrixFreeWorkspace{N,T<:Complex}
    work_real::Array{T,N}
    work_fourier::Array{T,N}
end

"""
    MatrixFreeWorkspace(resolution::NTuple{N, Int}, ::Type{T}=ComplexF64)

Create workspace arrays for the given resolution.
"""
function MatrixFreeWorkspace(
    resolution::NTuple{N,Int}, (::Type{T})=ComplexF64
) where {N,T<:Complex}
    MatrixFreeWorkspace{N,T}(zeros(T, resolution), zeros(T, resolution))
end

"""
    MatrixFreeWorkspace(ctx::FFTContext)

Create workspace arrays matching the FFT context.
"""
function MatrixFreeWorkspace(ctx::FFTContext{N,T,F,I}) where {N,T,F,I}
    MatrixFreeWorkspace(ctx.resolution, T)
end

# ============================================================================
# MatrixFreeOperator
# ============================================================================

"""
    MatrixFreeOperator{D<:Dimension, W<:WaveType, T<:Complex}

Matrix-free representation of the Hamiltonian operator H = RHS⁻¹ * LHS.

# Fields
- `solver`: Reference to the solver (contains geometry, basis, material_arrays)
- `k`: Current wave vector
- `ctx`: FFT context (shared, thread-safe for execution)
- `workspace`: Work arrays (thread-local)

# Thread Safety
The FFT plans in `ctx` can be safely shared across threads. However, each thread
must have its own `workspace` to avoid data races. For parallel k-point loops:

```julia
ctx = FFTContext(solver)  # Create once
Threads.@threads for k in k_points
    ws = MatrixFreeWorkspace(ctx)  # One per thread
    op = MatrixFreeOperator(solver, k, ctx, ws)
    # ... use op ...
end
```

# Grid Resolution Constraint

**Important**: The solver's grid resolution must be large enough to contain all
plane wave indices. The required resolution in each dimension is:

    resolution >= 2 * cutoff + 1

For example, with `cutoff=7` (default), the resolution must be at least 15 in each
dimension (e.g., 16×16 for 2D, 16×16×16 for 3D).

If the resolution is too small, plane wave coefficients outside the grid range
will be silently dropped during FFT operations, leading to incorrect results.
"""
struct MatrixFreeOperator{D<:Dimension,W<:WaveType,T<:Complex,N,F,I}
    solver::Solver{D,W}
    k::Vector{Float64}
    ctx::FFTContext{N,T,F,I}
    workspace::MatrixFreeWorkspace{N,T}
end

# Convenience accessors for backward compatibility
resolution(op::MatrixFreeOperator) = op.ctx.resolution
fft_plan(op::MatrixFreeOperator) = op.ctx.fft_plan
ifft_plan(op::MatrixFreeOperator) = op.ctx.ifft_plan
work_real(op::MatrixFreeOperator) = op.workspace.work_real
work_fourier(op::MatrixFreeOperator) = op.workspace.work_fourier

"""
    MatrixFreeOperator(solver::Solver, k, ctx::FFTContext, workspace::MatrixFreeWorkspace)

Create a matrix-free operator with explicit context and workspace (recommended for parallel use).
"""
function MatrixFreeOperator(
    solver::Solver{D,W},
    k::Vector{Float64},
    ctx::FFTContext{N,T,F,I},
    workspace::MatrixFreeWorkspace{N,T},
) where {D<:Dimension,W<:WaveType,N,T,F,I}
    MatrixFreeOperator{D,W,T,N,F,I}(solver, k, ctx, workspace)
end

function MatrixFreeOperator(
    solver::Solver{Dim1,W},
    k::Real,
    ctx::FFTContext{1,T,F,I},
    workspace::MatrixFreeWorkspace{1,T},
) where {W<:WaveType,T,F,I}
    MatrixFreeOperator{Dim1,W,T,1,F,I}(solver, [Float64(k)], ctx, workspace)
end

"""
    MatrixFreeOperator(solver::Solver, k)

Create a matrix-free operator (convenience method, creates new context and workspace).
For better performance in loops, create `FFTContext` once and reuse.
"""
function MatrixFreeOperator(solver::Solver{Dim2,W}, k::Vector{Float64}) where {W<:WaveType}
    ctx = FFTContext(solver.resolution, ComplexF64)
    workspace = MatrixFreeWorkspace(ctx)
    MatrixFreeOperator(solver, k, ctx, workspace)
end

function MatrixFreeOperator(solver::Solver{Dim1,W}, k::Real) where {W<:WaveType}
    ctx = FFTContext(solver.resolution, ComplexF64)
    workspace = MatrixFreeWorkspace(ctx)
    MatrixFreeOperator(solver, Float64(k), ctx, workspace)
end

# 1D convenience constructor accepting Vector{Float64} (for compatibility with solve interface)
function MatrixFreeOperator(solver::Solver{Dim1,W}, k::Vector{Float64}) where {W<:WaveType}
    MatrixFreeOperator(solver, k[1])
end

function MatrixFreeOperator(solver::Solver{Dim3,W}, k::Vector{Float64}) where {W<:WaveType}
    ctx = FFTContext(solver.resolution, ComplexF64)
    workspace = MatrixFreeWorkspace(ctx)
    MatrixFreeOperator(solver, k, ctx, workspace)
end

# ============================================================================
# Fourier space ↔ Grid space conversion
# ============================================================================

# Use same convention as convolution_matrix: fftshift(fft(fftshift(f)))
# For coefficients: place in shifted grid, then ifftshift(ifft(ifftshift))

"""
    fourier_to_grid!(grid, coeffs, basis, resolution)

Convert Fourier coefficients to real space grid.
Uses same convention as convolution_matrix.
"""
function fourier_to_grid!(
    grid::AbstractArray{T},
    coeffs::AbstractVector{T},
    basis::PlaneWaveBasis{Dim2},
    resolution::NTuple{2,Int},
) where {T}
    fill!(grid, zero(T))
    Nx, Ny = resolution
    cx = (Nx + 1) ÷ 2 + 1
    cy = (Ny + 1) ÷ 2 + 1

    # Place coefficients in shifted grid (centered at cx, cy)
    for (i, (p, q)) in enumerate(basis.indices)
        ix = cx + p
        iy = cy + q
        if 1 <= ix <= Nx && 1 <= iy <= Ny
            grid[ix, iy] = coeffs[i]
        end
    end

    # Transform to real space: ifftshift, then ifft, then fftshift
    # This is inverse of: fftshift(fft(fftshift(f))) / N
    grid .= fftshift(ifft(ifftshift(grid))) * (Nx * Ny)
    return grid
end

function fourier_to_grid!(
    grid::AbstractVector{T},
    coeffs::AbstractVector{T},
    basis::PlaneWaveBasis{Dim1},
    resolution::NTuple{1,Int},
) where {T}
    fill!(grid, zero(T))
    N = resolution[1]
    cx = (N + 1) ÷ 2 + 1

    for (i, (p,)) in enumerate(basis.indices)
        ix = cx + p
        if 1 <= ix <= N
            grid[ix] = coeffs[i]
        end
    end

    grid .= fftshift(ifft(ifftshift(grid))) * N
    return grid
end

function fourier_to_grid!(
    grid::AbstractArray{T,3},
    coeffs::AbstractVector{T},
    basis::PlaneWaveBasis{Dim3},
    resolution::NTuple{3,Int},
) where {T}
    fill!(grid, zero(T))
    Nx, Ny, Nz = resolution
    cx = (Nx + 1) ÷ 2 + 1
    cy = (Ny + 1) ÷ 2 + 1
    cz = (Nz + 1) ÷ 2 + 1

    for (i, (p, q, r)) in enumerate(basis.indices)
        ix = cx + p
        iy = cy + q
        iz = cz + r
        if 1 <= ix <= Nx && 1 <= iy <= Ny && 1 <= iz <= Nz
            grid[ix, iy, iz] = coeffs[i]
        end
    end

    grid .= fftshift(ifft(ifftshift(grid))) * (Nx * Ny * Nz)
    return grid
end

"""
    grid_to_fourier!(coeffs, grid, basis, resolution)

Convert real space grid to Fourier coefficients.
Uses same convention as convolution_matrix.
"""
function grid_to_fourier!(
    coeffs::AbstractVector{T},
    grid::AbstractArray{T},
    basis::PlaneWaveBasis{Dim2},
    resolution::NTuple{2,Int},
) where {T}
    Nx, Ny = resolution
    cx = (Nx + 1) ÷ 2 + 1
    cy = (Ny + 1) ÷ 2 + 1

    # FFT with same convention as convolution_matrix
    grid_fft = fftshift(fft(fftshift(grid))) / (Nx * Ny)

    for (i, (p, q)) in enumerate(basis.indices)
        ix = cx + p
        iy = cy + q
        if 1 <= ix <= Nx && 1 <= iy <= Ny
            coeffs[i] = grid_fft[ix, iy]
        else
            coeffs[i] = zero(T)
        end
    end

    return coeffs
end

function grid_to_fourier!(
    coeffs::AbstractVector{T},
    grid::AbstractVector{T},
    basis::PlaneWaveBasis{Dim1},
    resolution::NTuple{1,Int},
) where {T}
    N = resolution[1]
    cx = (N + 1) ÷ 2 + 1

    grid_fft = fftshift(fft(fftshift(grid))) / N

    for (i, (p,)) in enumerate(basis.indices)
        ix = cx + p
        if 1 <= ix <= N
            coeffs[i] = grid_fft[ix]
        else
            coeffs[i] = zero(T)
        end
    end

    return coeffs
end

function grid_to_fourier!(
    coeffs::AbstractVector{T},
    grid::AbstractArray{T,3},
    basis::PlaneWaveBasis{Dim3},
    resolution::NTuple{3,Int},
) where {T}
    Nx, Ny, Nz = resolution
    cx = (Nx + 1) ÷ 2 + 1
    cy = (Ny + 1) ÷ 2 + 1
    cz = (Nz + 1) ÷ 2 + 1

    grid_fft = fftshift(fft(fftshift(grid))) / (Nx * Ny * Nz)

    for (i, (p, q, r)) in enumerate(basis.indices)
        ix = cx + p
        iy = cy + q
        iz = cz + r
        if 1 <= ix <= Nx && 1 <= iy <= Ny && 1 <= iz <= Nz
            coeffs[i] = grid_fft[ix, iy, iz]
        else
            coeffs[i] = zero(T)
        end
    end

    return coeffs
end

# ============================================================================
# LHS operator application (H_LHS * v)
# ============================================================================

"""
    apply_lhs!(y, op::MatrixFreeOperator, x)

Apply LHS operator: y = LHS * x
For TE: LHS = Kx * ε⁻¹ * Kx + Ky * ε⁻¹ * Ky
"""
function apply_lhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim2,TEWave,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    k = op.k
    res = resolution(op)
    ε_inv = solver.material_arrays.ε_inv

    fill!(y, zero(T))
    N = basis.num_pw

    # Temporary arrays
    temp_grid = zeros(T, res)
    temp_coeffs = zeros(T, N)

    # For each direction (x, y)
    for dir in 1:2
        # Step 1: Apply K_dir in Fourier space (multiply by k_dir + G_dir)
        for i in eachindex(x)
            G = basis.G[i]
            K_dir = k[dir] + G[dir]
            temp_coeffs[i] = K_dir * x[i]
        end

        # Step 2: FFT to real space
        fourier_to_grid!(temp_grid, temp_coeffs, basis, res)

        # Step 3: Multiply by ε⁻¹ in real space
        temp_grid .*= ε_inv

        # Step 4: FFT back to Fourier space
        grid_to_fourier!(temp_coeffs, temp_grid, basis, res)

        # Step 5: Apply K_dir again in Fourier space
        for i in eachindex(y)
            G = basis.G[i]
            K_dir = k[dir] + G[dir]
            y[i] += K_dir * temp_coeffs[i]
        end
    end

    return y
end

"""
    apply_lhs!(y, op::MatrixFreeOperator{Dim1}, x)

Apply LHS operator for 1D: y = K * ε⁻¹ * K * x
"""
function apply_lhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim1,Photonic1D,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    k = op.k[1]
    res = resolution(op)
    ε_inv = solver.material_arrays.ε_inv

    N = basis.num_pw
    temp_grid = zeros(T, res)
    temp_coeffs = zeros(T, N)

    # Step 1: Apply K in Fourier space
    for i in eachindex(x)
        G = basis.G[i][1]
        K = k + G
        temp_coeffs[i] = K * x[i]
    end

    # Step 2: FFT to real space
    fourier_to_grid!(temp_grid, temp_coeffs, basis, res)

    # Step 3: Multiply by ε⁻¹
    temp_grid .*= ε_inv

    # Step 4: FFT back to Fourier space
    grid_to_fourier!(temp_coeffs, temp_grid, basis, res)

    # Step 5: Apply K again
    for i in eachindex(y)
        G = basis.G[i][1]
        K = k + G
        y[i] = K * temp_coeffs[i]
    end

    return y
end

# TM wave (photonic) - LHS = Kx * μ⁻¹ * Kx + Ky * μ⁻¹ * Ky
function apply_lhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim2,TMWave,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    k = op.k
    res = resolution(op)
    μ_inv = solver.material_arrays.μ_inv

    fill!(y, zero(T))
    N = basis.num_pw

    temp_grid = zeros(T, res)
    temp_coeffs = zeros(T, N)

    for dir in 1:2
        # Step 1: Apply K_dir in Fourier space
        for i in eachindex(x)
            G = basis.G[i]
            K_dir = k[dir] + G[dir]
            temp_coeffs[i] = K_dir * x[i]
        end

        # Step 2: FFT to real space
        fourier_to_grid!(temp_grid, temp_coeffs, basis, res)

        # Step 3: Multiply by μ⁻¹ in real space
        temp_grid .*= μ_inv

        # Step 4: FFT back to Fourier space
        grid_to_fourier!(temp_coeffs, temp_grid, basis, res)

        # Step 5: Apply K_dir again in Fourier space
        for i in eachindex(y)
            G = basis.G[i]
            K_dir = k[dir] + G[dir]
            y[i] += K_dir * temp_coeffs[i]
        end
    end

    return y
end

# SH wave (phononic)
function apply_lhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim2,SHWave,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    k = op.k
    res = resolution(op)
    C44 = solver.material_arrays.C44

    fill!(y, zero(T))
    N = basis.num_pw

    temp_grid = zeros(T, res)
    temp_coeffs = zeros(T, N)

    for dir in 1:2
        for i in eachindex(x)
            G = basis.G[i]
            K_dir = k[dir] + G[dir]
            temp_coeffs[i] = K_dir * x[i]
        end

        fourier_to_grid!(temp_grid, temp_coeffs, basis, res)
        temp_grid .*= C44
        grid_to_fourier!(temp_coeffs, temp_grid, basis, res)

        for i in eachindex(y)
            G = basis.G[i]
            K_dir = k[dir] + G[dir]
            y[i] += K_dir * temp_coeffs[i]
        end
    end

    return y
end

# Longitudinal 1D (phononic)
function apply_lhs!(
    y::AbstractVector{T},
    op::MatrixFreeOperator{Dim1,Longitudinal1D,T},
    x::AbstractVector{T},
) where {T}
    solver = op.solver
    basis = solver.basis
    k = op.k[1]
    res = resolution(op)
    C11 = solver.material_arrays.C11

    N = basis.num_pw
    temp_grid = zeros(T, res)
    temp_coeffs = zeros(T, N)

    for i in eachindex(x)
        G = basis.G[i][1]
        K = k + G
        temp_coeffs[i] = K * x[i]
    end

    fourier_to_grid!(temp_grid, temp_coeffs, basis, res)
    temp_grid .*= C11
    grid_to_fourier!(temp_coeffs, temp_grid, basis, res)

    for i in eachindex(y)
        G = basis.G[i][1]
        K = k + G
        y[i] = K * temp_coeffs[i]
    end

    return y
end

# PSV wave (phononic) - 2N components: [u_x; u_y]
# LHS = [K_xx K_xy; K_yx K_yy]
# K_xx = Kx * C11 * Kx + Ky * C44 * Ky
# K_yy = Kx * C44 * Kx + Ky * C11 * Ky
# K_xy = Kx * C12 * Ky + Ky * C44 * Kx
# K_yx = Kx * C44 * Ky + Ky * C12 * Kx
function apply_lhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim2,PSVWave,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    k = op.k
    res = resolution(op)
    C11 = solver.material_arrays.C11
    C12 = solver.material_arrays.C12
    C44 = solver.material_arrays.C44

    N = basis.num_pw
    x_x = @view x[1:N]      # u_x component
    x_y = @view x[(N + 1):2N]   # u_y component
    y_x = @view y[1:N]
    y_y = @view y[(N + 1):2N]

    fill!(y, zero(T))

    temp_grid = zeros(T, res)
    temp_coeffs = zeros(T, N)
    temp_coeffs2 = zeros(T, N)

    # y_x = K_xx * x_x + K_xy * x_y
    # K_xx = Kx * C11 * Kx + Ky * C44 * Ky

    # K_xx * x_x: Kx * C11 * Kx * x_x + Ky * C44 * Ky * x_x
    # Term 1: Kx * C11 * Kx * x_x
    for i in eachindex(x_x)
        G = basis.G[i]
        Kx = k[1] + G[1]
        temp_coeffs[i] = Kx * x_x[i]
    end
    fourier_to_grid!(temp_grid, temp_coeffs, basis, res)
    temp_grid .*= C11
    grid_to_fourier!(temp_coeffs, temp_grid, basis, res)
    for i in eachindex(y_x)
        G = basis.G[i]
        Kx = k[1] + G[1]
        y_x[i] += Kx * temp_coeffs[i]
    end

    # Term 2: Ky * C44 * Ky * x_x
    for i in eachindex(x_x)
        G = basis.G[i]
        Ky = k[2] + G[2]
        temp_coeffs[i] = Ky * x_x[i]
    end
    fourier_to_grid!(temp_grid, temp_coeffs, basis, res)
    temp_grid .*= C44
    grid_to_fourier!(temp_coeffs, temp_grid, basis, res)
    for i in eachindex(y_x)
        G = basis.G[i]
        Ky = k[2] + G[2]
        y_x[i] += Ky * temp_coeffs[i]
    end

    # K_xy * x_y: Kx * C12 * Ky * x_y + Ky * C44 * Kx * x_y
    # Term 3: Kx * C12 * Ky * x_y
    for i in eachindex(x_y)
        G = basis.G[i]
        Ky = k[2] + G[2]
        temp_coeffs[i] = Ky * x_y[i]
    end
    fourier_to_grid!(temp_grid, temp_coeffs, basis, res)
    temp_grid .*= C12
    grid_to_fourier!(temp_coeffs, temp_grid, basis, res)
    for i in eachindex(y_x)
        G = basis.G[i]
        Kx = k[1] + G[1]
        y_x[i] += Kx * temp_coeffs[i]
    end

    # Term 4: Ky * C44 * Kx * x_y
    for i in eachindex(x_y)
        G = basis.G[i]
        Kx = k[1] + G[1]
        temp_coeffs[i] = Kx * x_y[i]
    end
    fourier_to_grid!(temp_grid, temp_coeffs, basis, res)
    temp_grid .*= C44
    grid_to_fourier!(temp_coeffs, temp_grid, basis, res)
    for i in eachindex(y_x)
        G = basis.G[i]
        Ky = k[2] + G[2]
        y_x[i] += Ky * temp_coeffs[i]
    end

    # y_y = K_yx * x_x + K_yy * x_y
    # K_yx = Kx * C44 * Ky + Ky * C12 * Kx
    # K_yy = Kx * C44 * Kx + Ky * C11 * Ky

    # K_yx * x_x: Kx * C44 * Ky * x_x + Ky * C12 * Kx * x_x
    # Term 5: Kx * C44 * Ky * x_x
    for i in eachindex(x_x)
        G = basis.G[i]
        Ky = k[2] + G[2]
        temp_coeffs[i] = Ky * x_x[i]
    end
    fourier_to_grid!(temp_grid, temp_coeffs, basis, res)
    temp_grid .*= C44
    grid_to_fourier!(temp_coeffs, temp_grid, basis, res)
    for i in eachindex(y_y)
        G = basis.G[i]
        Kx = k[1] + G[1]
        y_y[i] += Kx * temp_coeffs[i]
    end

    # Term 6: Ky * C12 * Kx * x_x
    for i in eachindex(x_x)
        G = basis.G[i]
        Kx = k[1] + G[1]
        temp_coeffs[i] = Kx * x_x[i]
    end
    fourier_to_grid!(temp_grid, temp_coeffs, basis, res)
    temp_grid .*= C12
    grid_to_fourier!(temp_coeffs, temp_grid, basis, res)
    for i in eachindex(y_y)
        G = basis.G[i]
        Ky = k[2] + G[2]
        y_y[i] += Ky * temp_coeffs[i]
    end

    # K_yy * x_y: Kx * C44 * Kx * x_y + Ky * C11 * Ky * x_y
    # Term 7: Kx * C44 * Kx * x_y
    for i in eachindex(x_y)
        G = basis.G[i]
        Kx = k[1] + G[1]
        temp_coeffs[i] = Kx * x_y[i]
    end
    fourier_to_grid!(temp_grid, temp_coeffs, basis, res)
    temp_grid .*= C44
    grid_to_fourier!(temp_coeffs, temp_grid, basis, res)
    for i in eachindex(y_y)
        G = basis.G[i]
        Kx = k[1] + G[1]
        y_y[i] += Kx * temp_coeffs[i]
    end

    # Term 8: Ky * C11 * Ky * x_y
    for i in eachindex(x_y)
        G = basis.G[i]
        Ky = k[2] + G[2]
        temp_coeffs[i] = Ky * x_y[i]
    end
    fourier_to_grid!(temp_grid, temp_coeffs, basis, res)
    temp_grid .*= C11
    grid_to_fourier!(temp_coeffs, temp_grid, basis, res)
    for i in eachindex(y_y)
        G = basis.G[i]
        Ky = k[2] + G[2]
        y_y[i] += Ky * temp_coeffs[i]
    end

    return y
end

# ============================================================================
# RHS operator application
# ============================================================================

"""
    apply_rhs!(y, op::MatrixFreeOperator, x)

Apply RHS operator: y = RHS * x
"""
function apply_rhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim2,TEWave,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    res = resolution(op)
    μ = solver.material_arrays.μ

    N = basis.num_pw
    temp_grid = zeros(T, res)

    fourier_to_grid!(temp_grid, x, basis, res)
    temp_grid .*= μ
    grid_to_fourier!(y, temp_grid, basis, res)

    return y
end

# TM wave (photonic) - RHS = ε
function apply_rhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim2,TMWave,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    res = resolution(op)
    ε = solver.material_arrays.ε

    temp_grid = zeros(T, res)

    fourier_to_grid!(temp_grid, x, basis, res)
    temp_grid .*= ε
    grid_to_fourier!(y, temp_grid, basis, res)

    return y
end

function apply_rhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim2,SHWave,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    res = resolution(op)
    ρ = solver.material_arrays.ρ

    temp_grid = zeros(T, res)

    fourier_to_grid!(temp_grid, x, basis, res)
    temp_grid .*= ρ
    grid_to_fourier!(y, temp_grid, basis, res)

    return y
end

# PSV wave (phononic) - RHS = [ρ 0; 0 ρ] (block diagonal)
function apply_rhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim2,PSVWave,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    res = resolution(op)
    ρ = solver.material_arrays.ρ

    N = basis.num_pw
    x_x = @view x[1:N]
    x_y = @view x[(N + 1):2N]
    y_x = @view y[1:N]
    y_y = @view y[(N + 1):2N]

    temp_grid = zeros(T, res)

    # y_x = ρ * x_x
    fourier_to_grid!(temp_grid, x_x, basis, res)
    temp_grid .*= ρ
    grid_to_fourier!(y_x, temp_grid, basis, res)

    # y_y = ρ * x_y
    fourier_to_grid!(temp_grid, x_y, basis, res)
    temp_grid .*= ρ
    grid_to_fourier!(y_y, temp_grid, basis, res)

    return y
end

function apply_rhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim1,Photonic1D,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    res = resolution(op)
    μ = solver.material_arrays.μ

    temp_grid = zeros(T, res)

    fourier_to_grid!(temp_grid, x, basis, res)
    temp_grid .*= μ
    grid_to_fourier!(y, temp_grid, basis, res)

    return y
end

function apply_rhs!(
    y::AbstractVector{T},
    op::MatrixFreeOperator{Dim1,Longitudinal1D,T},
    x::AbstractVector{T},
) where {T}
    solver = op.solver
    basis = solver.basis
    res = resolution(op)
    ρ = solver.material_arrays.ρ

    temp_grid = zeros(T, res)

    fourier_to_grid!(temp_grid, x, basis, res)
    temp_grid .*= ρ
    grid_to_fourier!(y, temp_grid, basis, res)

    return y
end

# FullVectorEM (3D photonic) - LHS = curl × ε⁻¹ × curl
# H-field formulation:
# L_xx = Ky * ε⁻¹ * Ky + Kz * ε⁻¹ * Kz
# L_yy = Kx * ε⁻¹ * Kx + Kz * ε⁻¹ * Kz
# L_zz = Kx * ε⁻¹ * Kx + Ky * ε⁻¹ * Ky
# L_xy = -Ky * ε⁻¹ * Kx,  L_xz = -Kz * ε⁻¹ * Kx
# L_yx = -Kx * ε⁻¹ * Ky,  L_yz = -Kz * ε⁻¹ * Ky
# L_zx = -Kx * ε⁻¹ * Kz,  L_zy = -Ky * ε⁻¹ * Kz
function apply_lhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim3,FullVectorEM,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    k = op.k
    res = resolution(op)
    ε_inv = solver.material_arrays.ε_inv

    N = basis.num_pw
    x_x = @view x[1:N]
    x_y = @view x[(N + 1):2N]
    x_z = @view x[(2N + 1):3N]
    y_x = @view y[1:N]
    y_y = @view y[(N + 1):2N]
    y_z = @view y[(2N + 1):3N]

    fill!(y, zero(T))

    temp_grid = zeros(T, res)
    temp_coeffs = zeros(T, N)

    # Helper function: apply K_i * ε⁻¹ * K_j * input and accumulate to output with sign
    function apply_term!(output, input, i_dir, j_dir, sign)
        # Step 1: K_j * input
        for idx in eachindex(input)
            G = basis.G[idx]
            K_j = k[j_dir] + G[j_dir]
            temp_coeffs[idx] = K_j * input[idx]
        end
        # Step 2: FFT to real space
        fourier_to_grid!(temp_grid, temp_coeffs, basis, res)
        # Step 3: ε⁻¹ in real space
        temp_grid .*= ε_inv
        # Step 4: FFT back to Fourier
        grid_to_fourier!(temp_coeffs, temp_grid, basis, res)
        # Step 5: K_i and accumulate
        for idx in eachindex(output)
            G = basis.G[idx]
            K_i = k[i_dir] + G[i_dir]
            output[idx] += sign * K_i * temp_coeffs[idx]
        end
    end

    # y_x = L_xx * x_x + L_xy * x_y + L_xz * x_z
    # L_xx = Ky * ε⁻¹ * Ky + Kz * ε⁻¹ * Kz
    apply_term!(y_x, x_x, 2, 2, +1.0)  # Ky * ε⁻¹ * Ky * x_x
    apply_term!(y_x, x_x, 3, 3, +1.0)  # Kz * ε⁻¹ * Kz * x_x
    # L_xy = -Ky * ε⁻¹ * Kx
    apply_term!(y_x, x_y, 2, 1, -1.0)  # -Ky * ε⁻¹ * Kx * x_y
    # L_xz = -Kz * ε⁻¹ * Kx
    apply_term!(y_x, x_z, 3, 1, -1.0)  # -Kz * ε⁻¹ * Kx * x_z

    # y_y = L_yx * x_x + L_yy * x_y + L_yz * x_z
    # L_yx = -Kx * ε⁻¹ * Ky
    apply_term!(y_y, x_x, 1, 2, -1.0)  # -Kx * ε⁻¹ * Ky * x_x
    # L_yy = Kx * ε⁻¹ * Kx + Kz * ε⁻¹ * Kz
    apply_term!(y_y, x_y, 1, 1, +1.0)  # Kx * ε⁻¹ * Kx * x_y
    apply_term!(y_y, x_y, 3, 3, +1.0)  # Kz * ε⁻¹ * Kz * x_y
    # L_yz = -Kz * ε⁻¹ * Ky
    apply_term!(y_y, x_z, 3, 2, -1.0)  # -Kz * ε⁻¹ * Ky * x_z

    # y_z = L_zx * x_x + L_zy * x_y + L_zz * x_z
    # L_zx = -Kx * ε⁻¹ * Kz
    apply_term!(y_z, x_x, 1, 3, -1.0)  # -Kx * ε⁻¹ * Kz * x_x
    # L_zy = -Ky * ε⁻¹ * Kz
    apply_term!(y_z, x_y, 2, 3, -1.0)  # -Ky * ε⁻¹ * Kz * x_y
    # L_zz = Kx * ε⁻¹ * Kx + Ky * ε⁻¹ * Ky
    apply_term!(y_z, x_z, 1, 1, +1.0)  # Kx * ε⁻¹ * Kx * x_z
    apply_term!(y_z, x_z, 2, 2, +1.0)  # Ky * ε⁻¹ * Ky * x_z

    return y
end

# FullVectorEM (3D photonic) - RHS = diag(μ, μ, μ) (block diagonal)
# H-field formulation: RHS multiplies each component by μ independently
function apply_rhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim3,FullVectorEM,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    res = resolution(op)
    μ = solver.material_arrays.μ

    N = basis.num_pw
    x_x = @view x[1:N]
    x_y = @view x[(N + 1):2N]
    x_z = @view x[(2N + 1):3N]
    y_x = @view y[1:N]
    y_y = @view y[(N + 1):2N]
    y_z = @view y[(2N + 1):3N]

    temp_grid = zeros(T, res)

    # y_x = μ * x_x
    fourier_to_grid!(temp_grid, x_x, basis, res)
    temp_grid .*= μ
    grid_to_fourier!(y_x, temp_grid, basis, res)

    # y_y = μ * x_y
    fourier_to_grid!(temp_grid, x_y, basis, res)
    temp_grid .*= μ
    grid_to_fourier!(y_y, temp_grid, basis, res)

    # y_z = μ * x_z
    fourier_to_grid!(temp_grid, x_z, basis, res)
    temp_grid .*= μ
    grid_to_fourier!(y_z, temp_grid, basis, res)

    return y
end

# ============================================================================
# LinearMap interface
# ============================================================================

"""
    to_linear_map(op::MatrixFreeOperator)

Convert MatrixFreeOperator to a LinearMap for use with iterative solvers.
Returns LHS as a LinearMap.
"""
function to_linear_map_lhs(op::MatrixFreeOperator{D,W}) where {D,W}
    N = op.solver.basis.num_pw
    nc = ncomponents(op.solver.wave)
    dim = N * nc
    LinearMap{ComplexF64}(
        x -> begin
            y = zeros(ComplexF64, dim)
            apply_lhs!(y, op, x)
            y
        end, dim; ismutating=false, ishermitian=true
    )
end

"""
    to_linear_map_rhs(op::MatrixFreeOperator)

Convert RHS to a LinearMap.
"""
function to_linear_map_rhs(op::MatrixFreeOperator{D,W}) where {D,W}
    N = op.solver.basis.num_pw
    nc = ncomponents(op.solver.wave)
    dim = N * nc
    LinearMap{ComplexF64}(
        x -> begin
            y = zeros(ComplexF64, dim)
            apply_rhs!(y, op, x)
            y
        end, dim; ismutating=false, ishermitian=true
    )
end

# ============================================================================
# FullElastic (3D phononic) - Matrix-Free Implementation
# ============================================================================
# Formulation: -∇·(C:∇u) = ω² ρ u
# For isotropic material with Voigt notation (C11, C12, C44):
#
# LHS (stiffness matrix) 3x3 block structure:
# K_xx = Kx*C11*Kx + Ky*C44*Ky + Kz*C44*Kz  (3 terms)
# K_yy = Kx*C44*Kx + Ky*C11*Ky + Kz*C44*Kz  (3 terms)
# K_zz = Kx*C44*Kx + Ky*C44*Ky + Kz*C11*Kz  (3 terms)
# K_xy = Kx*C12*Ky + Ky*C44*Kx              (2 terms)
# K_xz = Kx*C12*Kz + Kz*C44*Kx              (2 terms)
# K_yx = Ky*C12*Kx + Kx*C44*Ky              (2 terms)
# K_yz = Ky*C12*Kz + Kz*C44*Ky              (2 terms)
# K_zx = Kz*C12*Kx + Kx*C44*Kz              (2 terms)
# K_zy = Kz*C12*Ky + Ky*C44*Kz              (2 terms)
# Total: 21 terms, 42 FFTs
#
# RHS: block diagonal ρ (3 FFT pairs)

# FullElastic (3D phononic) - RHS = diag(ρ, ρ, ρ)
function apply_rhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim3,FullElastic,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    res = resolution(op)
    ρ = solver.material_arrays.ρ

    N = basis.num_pw
    x_x = @view x[1:N]
    x_y = @view x[(N + 1):2N]
    x_z = @view x[(2N + 1):3N]
    y_x = @view y[1:N]
    y_y = @view y[(N + 1):2N]
    y_z = @view y[(2N + 1):3N]

    temp_grid = zeros(T, res)

    # y_x = ρ * x_x
    fourier_to_grid!(temp_grid, x_x, basis, res)
    temp_grid .*= ρ
    grid_to_fourier!(y_x, temp_grid, basis, res)

    # y_y = ρ * x_y
    fourier_to_grid!(temp_grid, x_y, basis, res)
    temp_grid .*= ρ
    grid_to_fourier!(y_y, temp_grid, basis, res)

    # y_z = ρ * x_z
    fourier_to_grid!(temp_grid, x_z, basis, res)
    temp_grid .*= ρ
    grid_to_fourier!(y_z, temp_grid, basis, res)

    return y
end

# FullElastic (3D phononic) - LHS = K (stiffness matrix)
function apply_lhs!(
    y::AbstractVector{T}, op::MatrixFreeOperator{Dim3,FullElastic,T}, x::AbstractVector{T}
) where {T}
    solver = op.solver
    basis = solver.basis
    k = op.k
    res = resolution(op)
    C11 = solver.material_arrays.C11
    C12 = solver.material_arrays.C12
    C44 = solver.material_arrays.C44

    N = basis.num_pw
    x_x = @view x[1:N]      # u_x component
    x_y = @view x[(N + 1):2N]   # u_y component
    x_z = @view x[(2N + 1):3N]  # u_z component
    y_x = @view y[1:N]
    y_y = @view y[(N + 1):2N]
    y_z = @view y[(2N + 1):3N]

    fill!(y, zero(T))

    temp_grid = zeros(T, res)
    temp_coeffs = zeros(T, N)

    # Helper function: apply K_i * Material * K_j * input and accumulate to output
    function apply_term!(output, input, i_dir, j_dir, material)
        # Step 1: K_j * input
        for idx in eachindex(input)
            G = basis.G[idx]
            K_j = k[j_dir] + G[j_dir]
            temp_coeffs[idx] = K_j * input[idx]
        end
        # Step 2: FFT to real space
        fourier_to_grid!(temp_grid, temp_coeffs, basis, res)
        # Step 3: Multiply by material in real space
        temp_grid .*= material
        # Step 4: FFT back to Fourier
        grid_to_fourier!(temp_coeffs, temp_grid, basis, res)
        # Step 5: K_i and accumulate
        for idx in eachindex(output)
            G = basis.G[idx]
            K_i = k[i_dir] + G[i_dir]
            output[idx] += K_i * temp_coeffs[idx]
        end
    end

    # ========================================================================
    # y_x = K_xx * x_x + K_xy * x_y + K_xz * x_z
    # ========================================================================
    # K_xx = Kx*C11*Kx + Ky*C44*Ky + Kz*C44*Kz (3 terms)
    apply_term!(y_x, x_x, 1, 1, C11)  # Kx * C11 * Kx * x_x
    apply_term!(y_x, x_x, 2, 2, C44)  # Ky * C44 * Ky * x_x
    apply_term!(y_x, x_x, 3, 3, C44)  # Kz * C44 * Kz * x_x

    # K_xy = Kx*C12*Ky + Ky*C44*Kx (2 terms)
    apply_term!(y_x, x_y, 1, 2, C12)  # Kx * C12 * Ky * x_y
    apply_term!(y_x, x_y, 2, 1, C44)  # Ky * C44 * Kx * x_y

    # K_xz = Kx*C12*Kz + Kz*C44*Kx (2 terms)
    apply_term!(y_x, x_z, 1, 3, C12)  # Kx * C12 * Kz * x_z
    apply_term!(y_x, x_z, 3, 1, C44)  # Kz * C44 * Kx * x_z

    # ========================================================================
    # y_y = K_yx * x_x + K_yy * x_y + K_yz * x_z
    # ========================================================================
    # K_yx = Ky*C12*Kx + Kx*C44*Ky (2 terms)
    apply_term!(y_y, x_x, 2, 1, C12)  # Ky * C12 * Kx * x_x
    apply_term!(y_y, x_x, 1, 2, C44)  # Kx * C44 * Ky * x_x

    # K_yy = Kx*C44*Kx + Ky*C11*Ky + Kz*C44*Kz (3 terms)
    apply_term!(y_y, x_y, 1, 1, C44)  # Kx * C44 * Kx * x_y
    apply_term!(y_y, x_y, 2, 2, C11)  # Ky * C11 * Ky * x_y
    apply_term!(y_y, x_y, 3, 3, C44)  # Kz * C44 * Kz * x_y

    # K_yz = Ky*C12*Kz + Kz*C44*Ky (2 terms)
    apply_term!(y_y, x_z, 2, 3, C12)  # Ky * C12 * Kz * x_z
    apply_term!(y_y, x_z, 3, 2, C44)  # Kz * C44 * Ky * x_z

    # ========================================================================
    # y_z = K_zx * x_x + K_zy * x_y + K_zz * x_z
    # ========================================================================
    # K_zx = Kz*C12*Kx + Kx*C44*Kz (2 terms)
    apply_term!(y_z, x_x, 3, 1, C12)  # Kz * C12 * Kx * x_x
    apply_term!(y_z, x_x, 1, 3, C44)  # Kx * C44 * Kz * x_x

    # K_zy = Kz*C12*Ky + Ky*C44*Kz (2 terms)
    apply_term!(y_z, x_y, 3, 2, C12)  # Kz * C12 * Ky * x_y
    apply_term!(y_z, x_y, 2, 3, C44)  # Ky * C44 * Kz * x_y

    # K_zz = Kx*C44*Kx + Ky*C44*Ky + Kz*C11*Kz (3 terms)
    apply_term!(y_z, x_z, 1, 1, C44)  # Kx * C44 * Kx * x_z
    apply_term!(y_z, x_z, 2, 2, C44)  # Ky * C44 * Ky * x_z
    apply_term!(y_z, x_z, 3, 3, C11)  # Kz * C11 * Kz * x_z

    return y
end
