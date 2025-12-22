# Last-Modified: 2025-12-15T21:30:00+09:00
# RSK Extension Tests
# These tests require ReducedShiftedKrylov.jl to be loaded

using Test
using PhoXonic
using ReducedShiftedKrylov  # This triggers the extension loading

@testset "RSK Extension Tests" begin

    # Helper function
    function setup_ldos_test()
        lat = square_lattice(1.0)
        air = Dielectric(1.0)
        rod = Dielectric(4.0)
        geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])
        solver = Solver(TEWave(), geo, (32, 32); cutoff=5)
        return solver
    end

    @testset "MatrixFreeGF compute_ldos" begin
        solver = setup_ldos_test()
        pos = [0.5, 0.5]
        ω_values = [0.3, 0.4, 0.5]
        k_points = [[0.0, 0.0]]

        # Direct method (reference)
        ldos_direct = compute_ldos(solver, pos, ω_values, k_points, DirectGF(); η=0.05)

        # MatrixFreeGF with ApproximateRHSInv (default)
        ldos_mf = compute_ldos(solver, pos, ω_values, k_points, MatrixFreeGF(); η=0.05)

        # MatrixFreeGF with CGRHSInv
        ldos_mf_cg = compute_ldos(
            solver, pos, ω_values, k_points, MatrixFreeGF(rhs_inv_method=CGRHSInv()); η=0.05
        )

        # All should find similar peak positions (even if absolute values differ)
        @test argmax(ldos_direct) == argmax(ldos_mf)
        @test argmax(ldos_direct) == argmax(ldos_mf_cg)
    end

    @testset "RSKGF compute_ldos" begin
        solver = setup_ldos_test()
        pos = [0.5, 0.5]
        ω_values = [0.3, 0.4, 0.5]
        k_points = [[0.0, 0.0]]

        # Direct method (reference)
        ldos_direct = compute_ldos(solver, pos, ω_values, k_points, DirectGF(); η=0.05)

        # RSKGF
        ldos_rsk = compute_ldos(solver, pos, ω_values, k_points, RSKGF(); η=0.05)

        # Peak positions should match
        @test argmax(ldos_direct) == argmax(ldos_rsk)
    end

    @testset "compute_greens_function with MatrixFreeGF" begin
        solver = setup_ldos_test()
        k = [0.1, 0.2]
        ω_values = [0.3, 0.4]
        N = solver.basis.num_pw
        source = randn(ComplexF64, N)

        # Direct
        G_direct = compute_greens_function(solver, k, ω_values, source, DirectGF(); η=0.05)

        # MatrixFreeGF
        G_mf = compute_greens_function(solver, k, ω_values, source, MatrixFreeGF(); η=0.05)

        @test length(G_direct) == length(ω_values)
        @test length(G_mf) == length(ω_values)
        @test all(length.(G_direct) .== N)
        @test all(length.(G_mf) .== N)
    end

    @testset "compute_greens_function with RSKGF" begin
        solver = setup_ldos_test()
        k = [0.1, 0.2]
        ω_values = [0.3, 0.4]
        N = solver.basis.num_pw
        source = randn(ComplexF64, N)

        # RSKGF
        G_rsk = compute_greens_function(solver, k, ω_values, source, RSKGF(); η=0.05)

        @test length(G_rsk) == length(ω_values)
        @test all(length.(G_rsk) .== N)
    end

    @testset "compute_dos with MatrixFreeGF" begin
        solver = setup_ldos_test()
        ω_values = [0.3, 0.4, 0.5]
        k_points = [[0.0, 0.0], [0.5, 0.0]]

        # Direct
        dos_direct = compute_dos(solver, ω_values, k_points, DirectGF(); η=0.05)

        # MatrixFreeGF (stochastic)
        dos_mf = compute_dos(solver, ω_values, k_points, MatrixFreeGF(); η=0.05, n_random=5)

        @test length(dos_direct) == length(ω_values)
        @test length(dos_mf) == length(ω_values)
        @test all(dos_direct .>= 0)
    end

    @testset "compute_dos with RSKGF" begin
        solver = setup_ldos_test()
        ω_values = [0.3, 0.4, 0.5]
        k_points = [[0.0, 0.0], [0.5, 0.0]]

        # RSKGF (stochastic)
        dos_rsk = compute_dos(solver, ω_values, k_points, RSKGF(); η=0.05, n_random=5)

        @test length(dos_rsk) == length(ω_values)
    end
end
