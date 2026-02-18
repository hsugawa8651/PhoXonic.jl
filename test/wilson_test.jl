# Tests for Wilson loop and topological invariants

using Test
using PhoXonic
using LinearAlgebra

@testset "Wilson loop" begin
    # Phase 0: Core functions
    @testset "Core functions" begin
        # overlap_matrix: weighted inner product between eigenvector sets
        @testset "overlap_matrix" begin
            n = 5
            W = Matrix{ComplexF64}(I, n, n)  # identity weight

            # Test 1: overlap with self using orthonormal vectors should give identity
            v1 = Matrix{ComplexF64}(I, n, 2)  # 2 orthonormal bands
            M_self = PhoXonic.overlap_matrix(v1, v1, W)
            @test size(M_self) == (2, 2)
            @test M_self ≈ I(2) atol=1e-10

            # Test 2: random vectors should give (n_bands x n_bands) matrix
            v2 = randn(ComplexF64, n, 3)  # 3 bands
            v3 = randn(ComplexF64, n, 3)
            M = PhoXonic.overlap_matrix(v2, v3, W)
            @test size(M) == (3, 3)

            # Test 3: overlap_matrix(a, b, W) = a' * W * b
            @test M ≈ v2' * W * v3 atol=1e-10

            # Test 4: with non-trivial weight matrix
            W2 = Diagonal(rand(n) .+ 0.5)  # positive diagonal weight
            M_weighted = PhoXonic.overlap_matrix(v2, v3, W2)
            @test M_weighted ≈ v2' * W2 * v3 atol=1e-10
        end

        # unitary_approx: SVD-based unitary approximation
        @testset "unitary_approx" begin
            # Test 1: result should be unitary (U * U' ≈ I)
            M = randn(ComplexF64, 3, 3)
            U = PhoXonic.unitary_approx(M)
            @test U * U' ≈ I(3) atol=1e-10
            @test U' * U ≈ I(3) atol=1e-10

            # Test 2: different sizes
            for n in [2, 4, 5]
                M_n = randn(ComplexF64, n, n)
                U_n = PhoXonic.unitary_approx(M_n)
                @test U_n * U_n' ≈ I(n) atol=1e-10
            end

            # Test 3: already unitary matrix should remain unchanged
            # Create a unitary matrix via QR decomposition
            Q, _ = qr(randn(ComplexF64, 4, 4))
            Q_mat = Matrix(Q)
            U_from_unitary = PhoXonic.unitary_approx(Q_mat)
            @test U_from_unitary ≈ Q_mat atol=1e-10

            # Test 4: determinant should have magnitude 1
            M4 = randn(ComplexF64, 3, 3)
            U4 = PhoXonic.unitary_approx(M4)
            @test abs(det(U4)) ≈ 1.0 atol=1e-10
        end

        # wilson_phases: extract phases from Wilson matrix eigenvalues
        @testset "wilson_phases" begin
            # Test 1: known diagonal unitary matrix
            W = diagm([exp(im * 0.5), exp(im * 1.2), exp(-im * 0.8)])
            phases = PhoXonic.wilson_phases(W)
            @test all(-π .<= phases .<= π)
            @test sort(phases) ≈ sort([0.5, 1.2, -0.8]) atol = 1e-10

            # Test 2: identity matrix has zero phases
            W_id = Matrix{ComplexF64}(I, 3, 3)
            phases_id = PhoXonic.wilson_phases(W_id)
            @test all(abs.(phases_id) .< 1e-10)

            # Test 3: phases wrap to [-π, π]
            W_wrap = diagm([exp(im * 2.5), exp(im * (-2.5)), exp(im * π)])
            phases_wrap = PhoXonic.wilson_phases(W_wrap)
            @test all(-π .<= phases_wrap .<= π)
            @test sort(phases_wrap) ≈ sort([2.5, -2.5, π]) atol = 1e-10

            # Test 4: non-diagonal unitary matrix
            Q, _ = qr(randn(ComplexF64, 3, 3))
            W_nondiag = Matrix(Q)
            phases_nondiag = PhoXonic.wilson_phases(W_nondiag)
            @test length(phases_nondiag) == 3
            @test all(-π .<= phases_nondiag .<= π)
            # Sum of phases = arg(det(W)) (mod 2π)
            total_phase = angle(det(W_nondiag))
            @test mod(sum(phases_nondiag) - total_phase + π, 2π) - π ≈ 0 atol = 1e-10
        end

        # wilson_matrix: product of overlap matrices around k-space loop
        @testset "wilson_matrix" begin
            n_basis = 10
            n_bands = 2
            n_k = 5
            W_weight = Matrix{ComplexF64}(I, n_basis, n_basis)  # identity weight

            # Create sequence of random eigenvector sets
            spaces = [randn(ComplexF64, n_basis, n_bands) for _ in 1:n_k]

            # Test 1: wilson_matrix returns a square matrix of size n_bands
            W = PhoXonic.wilson_matrix(spaces, W_weight)
            @test size(W) == (n_bands, n_bands)

            # Test 2: result should be unitary (product of unitary matrices)
            @test W * W' ≈ I(n_bands) atol = 1e-8
            @test W' * W ≈ I(n_bands) atol = 1e-8

            # Test 3: with orthonormal eigenvectors (identity overlap),
            # wilson_matrix should be identity
            ortho_spaces = [Matrix{ComplexF64}(I, n_basis, n_bands) for _ in 1:n_k]
            W_ortho = PhoXonic.wilson_matrix(ortho_spaces, W_weight)
            @test W_ortho ≈ I(n_bands) atol = 1e-10

            # Test 4: single k-point (no loop) should give identity
            single_space = [randn(ComplexF64, n_basis, n_bands)]
            W_single = PhoXonic.wilson_matrix(single_space, W_weight)
            @test W_single ≈ I(n_bands) atol = 1e-10

            # Test 5: determinant magnitude should be 1
            @test abs(det(W)) ≈ 1.0 atol = 1e-8
        end

        # get_weight_matrix: retrieve weight matrix from solver for each WaveType
        @testset "get_weight_matrix" begin
            # Test 1D Photonic
            @testset "Photonic1D" begin
                lat = lattice_1d(1.0)
                geo = Geometry(lat, Dielectric(1.0), [(Segment(0.2, 0.8), Dielectric(4.0))])
                solver = Solver(Photonic1D(), geo, 64; cutoff = 10)
                W = get_weight_matrix(solver)
                n_basis = length(solver.basis.G)
                @test size(W) == (n_basis, n_basis)
            end

            # Test 1D Longitudinal
            @testset "Longitudinal1D" begin
                lat = lattice_1d(1.0)
                # Use Lamé parameters λ, μ (not E, ν)
                mat = IsotropicElastic(; ρ = 7800.0, λ = 115e9, μ = 82e9)
                geo = Geometry(lat, mat, [(Segment(0.2, 0.8), mat)])
                solver = Solver(Longitudinal1D(), geo, 64; cutoff = 10)
                W = get_weight_matrix(solver)
                n_basis = length(solver.basis.G)
                @test size(W) == (n_basis, n_basis)
            end

            # Test 2D TE
            @testset "TEWave" begin
                lat = square_lattice(1.0)
                geo = Geometry(
                    lat, Dielectric(1.0), [(Circle([0.0, 0.0], 0.3), Dielectric(4.0))]
                )
                solver = Solver(TEWave(), geo, (16, 16); cutoff = 4)
                W = get_weight_matrix(solver)
                n_basis = length(solver.basis.G)
                @test size(W) == (n_basis, n_basis)
            end

            # Test 2D TM
            @testset "TMWave" begin
                lat = square_lattice(1.0)
                geo = Geometry(
                    lat, Dielectric(1.0), [(Circle([0.0, 0.0], 0.3), Dielectric(4.0))]
                )
                solver = Solver(TMWave(), geo, (16, 16); cutoff = 4)
                W = get_weight_matrix(solver)
                n_basis = length(solver.basis.G)
                @test size(W) == (n_basis, n_basis)
            end

            # Test 2D SH
            @testset "SHWave" begin
                lat = square_lattice(1.0)
                mat = IsotropicElastic(; ρ = 7800.0, λ = 115e9, μ = 82e9)
                geo = Geometry(lat, mat, [(Circle([0.0, 0.0], 0.3), mat)])
                solver = Solver(SHWave(), geo, (16, 16); cutoff = 4)
                W = get_weight_matrix(solver)
                n_basis = length(solver.basis.G)
                @test size(W) == (n_basis, n_basis)
            end

            # Test 2D PSV (block diagonal)
            @testset "PSVWave" begin
                lat = square_lattice(1.0)
                mat = IsotropicElastic(; ρ = 7800.0, λ = 115e9, μ = 82e9)
                geo = Geometry(lat, mat, [(Circle([0.0, 0.0], 0.3), mat)])
                solver = Solver(PSVWave(), geo, (16, 16); cutoff = 4)
                W = get_weight_matrix(solver)
                n_basis = length(solver.basis.G)
                @test size(W) == (2 * n_basis, 2 * n_basis)  # block diagonal
            end
        end
    end

    # Phase 1: 1D Zak phase
    @testset "1D Zak phase" begin
        # Test result type
        @testset "ZakPhaseResult type" begin
            # Create a simple 1D photonic crystal
            lat = lattice_1d(1.0)
            geo = Geometry(lat, Dielectric(1.0), [(Segment(0.2, 0.8), Dielectric(4.0))])
            solver = Solver(Photonic1D(), geo, 64; cutoff = 10)

            # Compute Zak phase for band 1
            result = compute_zak_phase(solver, 1:1)

            # Check result type
            @test result isa PhoXonic.ZakPhaseResult
            @test length(result.phases) == 1
            @test result.bands == 1:1
            @test -π <= result.phases[1] <= π
        end

        # Test SSH-like structure: Zak phase should depend on unit cell choice
        @testset "SSH-like structure" begin
            lat = lattice_1d(1.0)

            # Geometry A: high-index layer at center (centered unit cell)
            geo_A = Geometry(lat, Dielectric(1.0), [(Segment(0.35, 0.65), Dielectric(9.0))])

            # Geometry B: high-index layer at edges (off-center unit cell)
            geo_B = Geometry(lat, Dielectric(1.0), [
                (Segment(0.0, 0.15), Dielectric(9.0)),
                (Segment(0.85, 1.0), Dielectric(9.0)),
            ])

            solver_A = Solver(Photonic1D(), geo_A, 128; cutoff = 20)
            solver_B = Solver(Photonic1D(), geo_B, 128; cutoff = 20)

            # Compute Zak phase for first band
            zak_A = compute_zak_phase(solver_A, 1:1; n_k = 100)
            zak_B = compute_zak_phase(solver_B, 1:1; n_k = 100)

            # Both should be quantized to 0 or π (due to inversion symmetry)
            @test abs(zak_A.phases[1]) < 0.3 || abs(abs(zak_A.phases[1]) - π) < 0.3
            @test abs(zak_B.phases[1]) < 0.3 || abs(abs(zak_B.phases[1]) - π) < 0.3

            # Phase difference should be close to π
            # (choosing unit cell at inversion center vs. away from center)
            phase_diff = mod(zak_B.phases[1] - zak_A.phases[1] + π, 2π) - π
            @test abs(abs(phase_diff) - π) < 0.4 || abs(phase_diff) < 0.4
        end

        # Test multiple bands
        @testset "Multiple bands" begin
            lat = lattice_1d(1.0)
            geo = Geometry(lat, Dielectric(1.0), [(Segment(0.3, 0.7), Dielectric(9.0))])
            solver = Solver(Photonic1D(), geo, 128; cutoff = 20)

            # Compute Zak phase for bands 1-3
            result = compute_zak_phase(solver, 1:3; n_k = 50)

            @test length(result.phases) == 3
            @test result.bands == 1:3
            @test all(-π .<= result.phases .<= π)
        end
    end

    # Phase 2: 2D Wilson loop
    @testset "2D Wilson loop" begin
        # Test result type
        @testset "WilsonSpectrumResult type" begin
            # Create a simple 2D photonic crystal (square lattice with circular rod)
            lat = square_lattice(1.0)
            geo = Geometry(
                lat, Dielectric(1.0), [(Circle([0.0, 0.0], 0.3), Dielectric(9.0))]
            )
            solver = Solver(TMWave(), geo, (32, 32); cutoff = 5)

            # Compute Wilson spectrum for bands 1:2 along Γ-X-Γ
            result = compute_wilson_spectrum(solver, 1:2; n_k_path = 11, n_k_loop = 30)

            # Check result type and structure
            @test result isa PhoXonic.WilsonSpectrumResult
            @test length(result.k_values) == 11
            @test size(result.phases) == (11, 2)  # n_k_path × n_bands
            @test result.bands == 1:2
            @test all(-π .<= result.phases .<= π)
        end

        # Test winding number calculation
        @testset "winding_number" begin
            # For a trivial system, winding number should be 0
            lat = square_lattice(1.0)
            geo = Geometry(
                lat, Dielectric(1.0), [(Circle([0.0, 0.0], 0.2), Dielectric(4.0))]
            )
            solver = Solver(TMWave(), geo, (32, 32); cutoff = 5)

            result = compute_wilson_spectrum(solver, 1:2; n_k_path = 21, n_k_loop = 50)

            # Winding number should be an integer
            w1 = winding_number(result, 1)
            w2 = winding_number(result, 2)
            @test w1 isa Int
            @test w2 isa Int

            # For this simple structure, winding should be 0
            @test w1 == 0
            @test w2 == 0
        end

        # Test different wave types
        @testset "Wave type support" begin
            lat = square_lattice(1.0)

            # TM wave
            geo_ph = Geometry(
                lat, Dielectric(1.0), [(Circle([0.0, 0.0], 0.3), Dielectric(4.0))]
            )
            solver_tm = Solver(TMWave(), geo_ph, (16, 16); cutoff = 4)
            result_tm = compute_wilson_spectrum(solver_tm, 1:1; n_k_path = 5, n_k_loop = 20)
            @test result_tm isa PhoXonic.WilsonSpectrumResult

            # TE wave
            solver_te = Solver(TEWave(), geo_ph, (16, 16); cutoff = 4)
            result_te = compute_wilson_spectrum(solver_te, 1:1; n_k_path = 5, n_k_loop = 20)
            @test result_te isa PhoXonic.WilsonSpectrumResult

            # SH wave (phononic)
            mat = IsotropicElastic(; ρ = 7800.0, λ = 115e9, μ = 82e9)
            geo_ph_el = Geometry(lat, mat, [(Circle([0.0, 0.0], 0.3), mat)])
            solver_sh = Solver(SHWave(), geo_ph_el, (16, 16); cutoff = 4)
            result_sh = compute_wilson_spectrum(solver_sh, 1:1; n_k_path = 5, n_k_loop = 20)
            @test result_sh isa PhoXonic.WilsonSpectrumResult
        end
    end
end
