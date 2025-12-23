# Tests for TransverseEM wave type (3D photonic crystals with transverse projection)

using Test
using PhoXonic
using LinearAlgebra
using StaticArrays

@testset "TransverseEM" begin

    # ========================================================================
    # Phase 1: Polarization vector computation
    # ========================================================================
    @testset "compute_polarization_vectors" begin
        # Test case 1: K = [1, 0, 0] (along x-axis)
        @testset "K along x-axis" begin
            K = SVector(1.0, 0.0, 0.0)
            e1, e2 = PhoXonic.compute_polarization_vectors(K)

            # e1, e2 should be unit vectors
            @test norm(e1) ≈ 1.0
            @test norm(e2) ≈ 1.0

            # e1, e2 should be perpendicular to K
            @test abs(dot(e1, K)) < 1e-10
            @test abs(dot(e2, K)) < 1e-10

            # e1, e2 should be perpendicular to each other
            @test abs(dot(e1, e2)) < 1e-10
        end

        # Test case 2: K = [0, 0, 1] (along z-axis, special case)
        @testset "K along z-axis (special case)" begin
            K = SVector(0.0, 0.0, 1.0)
            e1, e2 = PhoXonic.compute_polarization_vectors(K)

            @test norm(e1) ≈ 1.0
            @test norm(e2) ≈ 1.0
            @test abs(dot(e1, K)) < 1e-10
            @test abs(dot(e2, K)) < 1e-10
            @test abs(dot(e1, e2)) < 1e-10
        end

        # Test case 3: K = [1, 1, 1] (general case)
        @testset "K general direction" begin
            K = SVector(1.0, 1.0, 1.0)
            e1, e2 = PhoXonic.compute_polarization_vectors(K)

            @test norm(e1) ≈ 1.0
            @test norm(e2) ≈ 1.0
            @test abs(dot(e1, K)) < 1e-10
            @test abs(dot(e2, K)) < 1e-10
            @test abs(dot(e1, e2)) < 1e-10
        end

        # Test case 4: K ≈ 0 (Γ point)
        @testset "K near zero (Gamma point)" begin
            K = SVector(1e-12, 1e-12, 1e-12)
            e1, e2 = PhoXonic.compute_polarization_vectors(K)

            # Should return some valid orthonormal pair
            @test norm(e1) ≈ 1.0
            @test norm(e2) ≈ 1.0
            @test abs(dot(e1, e2)) < 1e-10
        end

        # Test case 5: Various random directions
        @testset "Random directions" begin
            for _ in 1:10
                K = SVector(randn(), randn(), randn())
                if norm(K) > 1e-10
                    e1, e2 = PhoXonic.compute_polarization_vectors(K)

                    @test norm(e1) ≈ 1.0 atol=1e-10
                    @test norm(e2) ≈ 1.0 atol=1e-10
                    @test abs(dot(e1, K / norm(K))) < 1e-10
                    @test abs(dot(e2, K / norm(K))) < 1e-10
                    @test abs(dot(e1, e2)) < 1e-10
                end
            end
        end
    end

    # ========================================================================
    # Phase 2: Polarization basis matrix
    # ========================================================================
    @testset "build_polarization_basis" begin
        # Simple FCC lattice setup
        a = 1.0
        lat = fcc_lattice(a)
        air = Dielectric(1.0)
        dielectric = Dielectric(11.56)

        geo = Geometry(lat, air, [(Sphere([0.0, 0.0, 0.0], 0.2), dielectric),])

        resolution = (4, 4, 4)
        cutoff = 2

        # Need FullVectorEM solver to get basis (TransverseEM not yet implemented)
        solver = Solver(FullVectorEM(), geo, resolution, DenseMethod(); cutoff=cutoff)
        N = solver.basis.num_pw

        k = [0.5, 0.3, 0.1]  # General k-point

        @testset "Matrix dimensions" begin
            P = PhoXonic.build_polarization_basis(solver.basis, k)
            @test size(P) == (3N, 2N)
        end

        @testset "Columns are orthonormal" begin
            P = PhoXonic.build_polarization_basis(solver.basis, k)
            # P'P should be identity (2N × 2N)
            PtP = P' * P
            @test PtP ≈ I(2N) atol=1e-10
        end
    end

    # ========================================================================
    # Phase 3: Matrix construction
    # ========================================================================
    @testset "build_matrices for TransverseEM" begin
        # Setup: FCC lattice with dielectric spheres
        a = 1.0
        lat = fcc_lattice(a)
        air = Dielectric(1.0)
        dielectric = Dielectric(11.56)

        geo = Geometry(lat, air, [(Sphere([0.0, 0.0, 0.0], 0.2), dielectric),])

        resolution = (4, 4, 4)
        cutoff = 2

        solver = Solver(TransverseEM(), geo, resolution, DenseMethod(); cutoff=cutoff)
        N = solver.basis.num_pw

        k = [0.5, 0.3, 0.1]  # General k-point

        @testset "Matrix dimensions" begin
            LHS, RHS = build_matrices(solver, k)
            @test size(LHS) == (2N, 2N)
            @test size(RHS) == (2N, 2N)
        end

        @testset "Matrices are Hermitian" begin
            LHS, RHS = build_matrices(solver, k)
            # Check Hermiticity: A ≈ A' up to numerical precision
            # Use relative tolerance for large matrices
            @test norm(LHS - LHS') / norm(LHS) < 1e-10
            @test norm(RHS - RHS') / norm(RHS) < 1e-10
        end

        @testset "RHS is positive definite" begin
            LHS, RHS = build_matrices(solver, k)
            # Check that all eigenvalues of RHS are positive
            eigvals_rhs = eigvals(Hermitian(RHS))
            @test all(eigvals_rhs .> 0)
        end

        @testset "Uniform medium analytical solution" begin
            # Uniform dielectric: ε = 4.0 (n = 2)
            # Expected: ω = c|k|/n = |k|/2
            ε_uniform = 4.0
            air_uniform = Dielectric(ε_uniform)
            geo_uniform = Geometry(lat, air_uniform)  # No inclusions

            solver_uniform = Solver(
                TransverseEM(), geo_uniform, resolution, DenseMethod(); cutoff=cutoff
            )

            # Use a simple k-point
            k_test = [0.5, 0.0, 0.0]
            k_mag = norm(k_test)
            n = sqrt(ε_uniform)  # refractive index
            ω_expected = k_mag / n  # ω = c|k|/n, c=1 in normalized units

            freqs, _ = solve(solver_uniform, k_test; bands=1:4)

            # First two modes should be doubly degenerate transverse modes
            # with ω = |k|/n
            @test freqs[1] ≈ ω_expected atol=0.01
            @test freqs[2] ≈ ω_expected atol=0.01
        end
    end

    # ========================================================================
    # Phase 4: Solver integration
    # ========================================================================
    @testset "solve for TransverseEM" begin
        # Setup: FCC lattice with dielectric spheres
        a = 1.0
        lat = fcc_lattice(a)
        air = Dielectric(1.0)
        dielectric = Dielectric(11.56)

        geo = Geometry(lat, air, [(Sphere([0.0, 0.0, 0.0], 0.2), dielectric),])

        resolution = (4, 4, 4)
        cutoff = 2
        solver = Solver(TransverseEM(), geo, resolution, DenseMethod(); cutoff=cutoff)
        N = solver.basis.num_pw

        @testset "Return types and sizes" begin
            k = [0.5, 0.3, 0.1]
            freqs, vecs = solve(solver, k; bands=1:6)

            # Frequencies should be Vector{Float64}
            @test freqs isa Vector{Float64}
            @test length(freqs) == 6

            # Eigenvectors should be Matrix{ComplexF64} with size (2N, nbands)
            @test vecs isa Matrix{ComplexF64}
            @test size(vecs) == (2N, 6)

            # Frequencies should be non-negative and sorted
            @test all(freqs .>= 0)
            @test issorted(freqs)
        end

        @testset "No spurious longitudinal modes" begin
            # At a general k-point, all modes should have ω > 0
            k = [0.3, 0.2, 0.1]
            freqs, _ = solve(solver, k; bands=1:10)

            # All frequencies should be significantly positive
            # (spurious modes in FullVectorEM would have ω ≈ 0)
            @test all(freqs .> 0.01)
        end

        @testset "Comparison with FullVectorEM transverse modes" begin
            # For a structure, TransverseEM should give the same physical modes
            # as FullVectorEM, but without the spurious longitudinal modes
            k = [0.3, 0.2, 0.1]

            # TransverseEM
            solver_trans = Solver(
                TransverseEM(), geo, resolution, DenseMethod(); cutoff=cutoff
            )
            freqs_trans, _ = solve(solver_trans, k; bands=1:4)

            # FullVectorEM (filter out spurious modes with ω ≈ 0)
            solver_full = Solver(
                FullVectorEM(), geo, resolution, DenseMethod(); cutoff=cutoff
            )
            freqs_full, _ = solve(solver_full, k; bands=1:10)
            freqs_full_physical = filter(f -> f > 0.01, freqs_full)

            # First few physical modes should match
            for i in 1:min(4, length(freqs_full_physical))
                @test freqs_trans[i] ≈ freqs_full_physical[i] rtol=0.01
            end
        end
    end

    # ========================================================================
    # Phase 5: MPB validation
    # ========================================================================
    @testset "FCC photonic crystal" begin
        # FCC lattice with ε=11.56 spheres in air, r=0.25a
        a = 1.0
        lat = fcc_lattice(a)
        air = Dielectric(1.0)
        dielectric = Dielectric(11.56)  # Silicon-like

        geo = Geometry(lat, air, [(Sphere([0.0, 0.0, 0.0], 0.25), dielectric),])

        # Use moderate resolution for tests
        resolution = (8, 8, 8)
        cutoff = 4

        solver = Solver(TransverseEM(), geo, resolution, DenseMethod(); cutoff=cutoff)

        @testset "Physical frequencies at high-symmetry points" begin
            # Test at a k-point along Γ-X direction (using Cartesian coordinates)
            # k = [π, 0, 0] (halfway to zone boundary)
            k_test = [π, 0.0, 0.0]
            freqs, _ = solve(solver, k_test; bands=1:6)

            # All eigenvalues should be positive (physical frequencies)
            @test all(freqs .> 0)

            # Normalized frequencies should be in reasonable range (0.1-1.0)
            freqs_normalized = freqs ./ (2π)
            @test all(freqs_normalized .> 0.1)
            @test all(freqs_normalized .< 2.0)
        end

        @testset "Band structure computation" begin
            # Compute bands along a k-path
            kpath = simple_kpath_fcc(; npoints=5)
            bands_result = compute_bands(solver, kpath; bands=1:4)

            # Verify output structure
            @test nbands(bands_result) == 4
            @test nkpoints(bands_result) > 0

            # All frequencies should be non-negative
            all_freqs = frequencies(bands_result)
            @test all(f -> all(f .>= 0), all_freqs)
        end
    end

    @testset "Resolution convergence" begin
        # Test that higher resolution gives more accurate results
        a = 1.0
        lat = fcc_lattice(a)
        air = Dielectric(1.0)
        dielectric = Dielectric(11.56)

        geo = Geometry(lat, air, [(Sphere([0.0, 0.0, 0.0], 0.25), dielectric),])

        k_test = [0.5, 0.25, 0.0]  # General k-point

        # Test at different resolutions
        resolutions = [(4, 4, 4), (6, 6, 6), (8, 8, 8)]
        cutoff = 3

        freqs_at_res = Float64[]
        for res in resolutions
            solver = Solver(TransverseEM(), geo, res, DenseMethod(); cutoff=cutoff)
            freqs, _ = solve(solver, k_test; bands=1:2)
            push!(freqs_at_res, freqs[1])
        end

        # Frequencies should converge (differences should decrease)
        diff1 = abs(freqs_at_res[2] - freqs_at_res[1])
        diff2 = abs(freqs_at_res[3] - freqs_at_res[2])
        @test diff2 < diff1  # Convergence: second difference should be smaller
    end
end
