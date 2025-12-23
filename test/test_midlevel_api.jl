# Test Mid-level API
# TDD: These tests should FAIL initially, then PASS after implementation

using Test
using PhoXonic
using LinearAlgebra

@testset "Mid-level API" begin
    # Setup: Create a simple 2D photonic crystal solver
    lat = square_lattice(1.0)
    air = Dielectric(1.0)
    rod = Dielectric(8.9)
    geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])

    solver_tm = Solver(TMWave(), geo, (16, 16); cutoff=5)
    solver_te = Solver(TEWave(), geo, (16, 16); cutoff=5)

    k = [0.1, 0.2]

    @testset "AbstractSolver type hierarchy" begin
        @test Solver <: AbstractSolver
        @test solver_tm isa AbstractSolver
        @test solver_te isa AbstractSolver
    end

    @testset "solve_at_k_with_vectors" begin
        # Test return types
        ω, vecs = solve_at_k_with_vectors(solver_tm, k, DenseMethod(); bands=1:4)

        @test ω isa Vector{Float64}
        @test vecs isa Matrix{ComplexF64}
        @test length(ω) == 4
        @test size(vecs, 2) == 4

        # Matrix dimension should match
        N = matrix_dimension(solver_tm)
        @test size(vecs, 1) == N

        # Frequencies should be non-negative and sorted
        @test all(ω .>= 0)
        @test issorted(ω)

        # Test with different methods
        ω_krylov, vecs_krylov = solve_at_k_with_vectors(
            solver_tm, k, KrylovKitMethod(); bands=1:4
        )
        @test ω_krylov isa Vector{Float64}
        @test vecs_krylov isa Matrix{ComplexF64}

        ω_lobpcg, vecs_lobpcg = solve_at_k_with_vectors(
            solver_tm, k, LOBPCGMethod(); bands=1:4
        )
        @test ω_lobpcg isa Vector{Float64}
        @test vecs_lobpcg isa Matrix{ComplexF64}
    end

    @testset "build_matrices" begin
        LHS, RHS = build_matrices(solver_tm, k)

        @test LHS isa Matrix{ComplexF64}
        @test RHS isa Matrix{ComplexF64}
        @test size(LHS) == size(RHS)

        N = matrix_dimension(solver_tm)
        @test size(LHS) == (N, N)

        # Test TE as well
        LHS_te, RHS_te = build_matrices(solver_te, k)
        @test LHS_te isa Matrix{ComplexF64}
        @test RHS_te isa Matrix{ComplexF64}
    end

    @testset "get_weight_matrix" begin
        W = get_weight_matrix(solver_tm)

        @test W isa Matrix{ComplexF64}
        N = matrix_dimension(solver_tm)
        @test size(W) == (N, N)

        # W should equal RHS from build_matrices
        _, RHS = build_matrices(solver_tm, k)
        @test W ≈ RHS atol=1e-6

        # Test TE
        W_te = get_weight_matrix(solver_te)
        @test W_te isa Matrix{ComplexF64}
    end

    @testset "Eigenvector orthonormality with W" begin
        ω, vecs = solve_at_k_with_vectors(solver_tm, k, DenseMethod(); bands=1:4)
        W = get_weight_matrix(solver_tm)

        # vecs' * W * vecs should be identity
        overlap = vecs' * W * vecs
        @test overlap ≈ I(4) atol=1e-5

        # Test with TE
        ω_te, vecs_te = solve_at_k_with_vectors(solver_te, k, DenseMethod(); bands=1:4)
        W_te = get_weight_matrix(solver_te)
        overlap_te = vecs_te' * W_te * vecs_te
        @test overlap_te ≈ I(4) atol=1e-5
    end

    @testset "Error fallbacks" begin
        # Invalid solver type
        @test_throws ErrorException solve_at_k_with_vectors(123, k, DenseMethod())
        @test_throws ErrorException solve_at_k_with_vectors("invalid", k, DenseMethod())

        # Invalid k type
        @test_throws ErrorException solve_at_k_with_vectors(
            solver_tm, "invalid", DenseMethod()
        )

        # Invalid method type
        @test_throws ErrorException solve_at_k_with_vectors(solver_tm, k, 123)

        # build_matrices errors
        @test_throws ErrorException build_matrices(nothing, k)
        @test_throws ErrorException build_matrices(123, k)

        # get_weight_matrix errors
        @test_throws ErrorException get_weight_matrix(:invalid)
        @test_throws ErrorException get_weight_matrix(nothing)
    end

    @testset "2D Phononic waves" begin
        # Setup phononic solver
        epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
        steel = IsotropicElastic(ρ=7800.0, λ=1.15e11, μ=8.28e10)
        geo_ph = Geometry(lat, epoxy, [(Circle([0.0, 0.0], 0.3), steel)])

        solver_sh = Solver(SHWave(), geo_ph, (16, 16); cutoff=5)
        solver_psv = Solver(PSVWave(), geo_ph, (16, 16); cutoff=5)

        # SH wave (higher tolerance for high-contrast materials)
        ω_sh, vecs_sh = solve_at_k_with_vectors(solver_sh, k, DenseMethod(); bands=1:4)
        W_sh = get_weight_matrix(solver_sh)
        @test vecs_sh' * W_sh * vecs_sh ≈ I(4) atol=1e-3

        # PSV wave (higher tolerance for high-contrast materials)
        ω_psv, vecs_psv = solve_at_k_with_vectors(solver_psv, k, DenseMethod(); bands=1:4)
        W_psv = get_weight_matrix(solver_psv)
        @test vecs_psv' * W_psv * vecs_psv ≈ I(4) atol=1e-3
    end

    @testset "1D waves" begin
        # 1D Photonic
        lat_1d = lattice_1d(1.0)
        mat1 = Dielectric(1.0)
        mat2 = Dielectric(4.0)
        geo_1d = Geometry(lat_1d, mat1, [(Segment(0.0, 0.3), mat2)])

        solver_1d_ph = Solver(Photonic1D(), geo_1d, 32; cutoff=10)

        k_1d = 0.3
        ω_1d, vecs_1d = solve_at_k_with_vectors(
            solver_1d_ph, k_1d, DenseMethod(); bands=1:4
        )
        W_1d = get_weight_matrix(solver_1d_ph)
        @test vecs_1d' * W_1d * vecs_1d ≈ I(4) atol=1e-6

        # 1D Elastic
        al = IsotropicElastic(ρ=2700.0, λ=5.1e10, μ=2.6e10)
        steel_1d = IsotropicElastic(ρ=7800.0, λ=1.15e11, μ=8.28e10)
        geo_1d_el = Geometry(lat_1d, al, [(Segment(0.0, 0.5), steel_1d)])

        solver_1d_el = Solver(Longitudinal1D(), geo_1d_el, 32; cutoff=10)

        ω_1d_el, vecs_1d_el = solve_at_k_with_vectors(
            solver_1d_el, k_1d, DenseMethod(); bands=1:4
        )
        W_1d_el = get_weight_matrix(solver_1d_el)
        @test vecs_1d_el' * W_1d_el * vecs_1d_el ≈ I(4) atol=1e-5
    end
end
