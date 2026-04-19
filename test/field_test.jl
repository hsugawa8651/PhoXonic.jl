# test/field_test.jl
# Field reconstruction tests for PhoXonic.jl

using Test
using PhoXonic
using LinearAlgebra

@testset "Field Reconstruction" begin

    # =========================================================================
    # Phase 1A-1: 1D scalar field
    # =========================================================================
    @testset "1D Photonic" begin
        # Setup: 1D photonic crystal
        lat = lattice_1d(1.0)
        geo = Geometry(lat, Dielectric(1.0), [
            (Segment(0.0, 0.3), Dielectric(9.0))
        ])
        solver = Solver(Photonic1D(), geo, (64,); cutoff=5)

        # Solve at k = 0
        freqs, vecs = solve_at_k_with_vectors(solver, [0.0], DenseMethod())

        @testset "reconstruct_field returns correct type and size" begin
            field = reconstruct_field(solver, vecs[:, 1])

            @test field isa Vector{ComplexF64}
            @test length(field) == 64
        end

        @testset "reconstruct_field returns non-zero field" begin
            field = reconstruct_field(solver, vecs[:, 1])

            @test !all(iszero, field)
            @test maximum(abs.(field)) > 0
        end

        @testset "custom grid size" begin
            field = reconstruct_field(solver, vecs[:, 1]; grid=(128,))

            @test length(field) == 128
        end

        @testset "multiple bands have different fields" begin
            field1 = reconstruct_field(solver, vecs[:, 1])
            field2 = reconstruct_field(solver, vecs[:, 2])

            # Fields should be different (unless they are the same band)
            @test !isapprox(field1, field2; rtol=0.1)
        end
    end

    @testset "1D Longitudinal (Phononic)" begin
        lat = lattice_1d(1.0)
        geo = Geometry(lat, IsotropicElastic(ρ=1000.0, λ=1e9, μ=1e9), [
            (Segment(0.0, 0.3), IsotropicElastic(ρ=2000.0, λ=2e9, μ=2e9))
        ])
        solver = Solver(Longitudinal1D(), geo, (64,); cutoff=5)
        freqs, vecs = solve_at_k_with_vectors(solver, [0.0], DenseMethod())

        @testset "basic reconstruction" begin
            field = reconstruct_field(solver, vecs[:, 1])

            @test field isa Vector{ComplexF64}
            @test length(field) == 64
            @test !all(iszero, field)
        end
    end

    # =========================================================================
    # Phase 1A-2: 2D scalar field
    # =========================================================================
    @testset "2D TM (scalar)" begin
        lat = square_lattice(1.0)
        geo = Geometry(lat, Dielectric(1.0), [
            (Circle([0.5, 0.5], 0.2), Dielectric(9.0))
        ])
        solver = Solver(TMWave(), geo, (32, 32); cutoff=5)
        freqs, vecs = solve_at_k_with_vectors(solver, [0.0, 0.0], DenseMethod())

        @testset "reconstruct_field returns correct type and size" begin
            field = reconstruct_field(solver, vecs[:, 1])

            @test field isa Matrix{ComplexF64}
            @test size(field) == (32, 32)
        end

        @testset "non-zero field" begin
            field = reconstruct_field(solver, vecs[:, 1])

            @test !all(iszero, field)
        end

        @testset "custom grid size" begin
            field = reconstruct_field(solver, vecs[:, 1]; grid=(64, 64))

            @test size(field) == (64, 64)
        end
    end

    @testset "2D SH (scalar phononic)" begin
        lat = square_lattice(1.0)
        geo = Geometry(lat, IsotropicElastic(ρ=1000.0, λ=1e9, μ=1e9), [
            (Circle([0.5, 0.5], 0.2), IsotropicElastic(ρ=2000.0, λ=2e9, μ=2e9))
        ])
        solver = Solver(SHWave(), geo, (32, 32); cutoff=5)
        freqs, vecs = solve_at_k_with_vectors(solver, [0.0, 0.0], DenseMethod())

        @testset "basic reconstruction" begin
            field = reconstruct_field(solver, vecs[:, 1])

            @test field isa Matrix{ComplexF64}
            @test size(field) == (32, 32)
        end
    end

    # =========================================================================
    # Phase 1A-3: 2D vector field
    # =========================================================================
    @testset "2D TE (scalar: Hz)" begin
        lat = square_lattice(1.0)
        geo = Geometry(lat, Dielectric(1.0), [
            (Circle([0.5, 0.5], 0.2), Dielectric(9.0))
        ])
        solver = Solver(TEWave(), geo, (32, 32); cutoff=5)
        freqs, vecs = solve_at_k_with_vectors(solver, [0.1, 0.1], DenseMethod())  # non-zero k

        @testset "reconstruct_field returns scalar for TE" begin
            field = reconstruct_field(solver, vecs[:, 1])

            # TEWave has ncomponents=1 so returns scalar
            @test field isa Matrix{ComplexF64}
            @test size(field) == (32, 32)
        end
    end

    @testset "2D PSV (vector: ux, uy)" begin
        lat = square_lattice(1.0)
        geo = Geometry(lat, IsotropicElastic(ρ=1000.0, λ=1e9, μ=1e9), [
            (Circle([0.5, 0.5], 0.2), IsotropicElastic(ρ=2000.0, λ=2e9, μ=2e9))
        ])
        solver = Solver(PSVWave(), geo, (32, 32); cutoff=5)
        freqs, vecs = solve_at_k_with_vectors(solver, [0.1, 0.1], DenseMethod())

        @testset "reconstruct_field returns tuple for 2-component field" begin
            field = reconstruct_field(solver, vecs[:, 1])

            # PSVWave has ncomponents=2 so returns 2-element tuple
            @test field isa Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}
            @test length(field) == 2
            @test size(field[1]) == (32, 32)
            @test size(field[2]) == (32, 32)
        end

        @testset "both components are non-zero (at non-zero k)" begin
            field = reconstruct_field(solver, vecs[:, 1])

            @test !all(iszero, field[1])
            @test !all(iszero, field[2])
        end
    end

    # =========================================================================
    # Phase 1A-4: 3D field
    # =========================================================================
    @testset "3D TransverseEM" begin
        lat = cubic_lattice(1.0)
        geo = Geometry(lat, Dielectric(1.0), [
            (Sphere([0.5, 0.5, 0.5], 0.2), Dielectric(9.0))
        ])
        solver = Solver(TransverseEM(), geo, (16, 16, 16); cutoff=3)
        freqs, vecs = solve_at_k_with_vectors(solver, [0.1, 0.1, 0.1], DenseMethod())

        @testset "reconstruct_field returns tuple for 2-component 3D field" begin
            field = reconstruct_field(solver, vecs[:, 1])

            # TransverseEM has ncomponents=2
            @test field isa Tuple{Array{ComplexF64,3}, Array{ComplexF64,3}}
            @test length(field) == 2
            @test size(field[1]) == (16, 16, 16)
            @test size(field[2]) == (16, 16, 16)
        end
    end

    # =========================================================================
    # Phase 1A-5: get_epsilon_field / get_material_field
    # =========================================================================
    @testset "get_epsilon_field" begin
        @testset "2D Photonic" begin
            lat = square_lattice(1.0)
            geo = Geometry(lat, Dielectric(1.0), [
                (Circle([0.5, 0.5], 0.2), Dielectric(9.0))
            ])
            solver = Solver(TMWave(), geo, (32, 32); cutoff=5)

            eps = get_epsilon_field(solver)

            @test eps isa Matrix{Float64}
            @test size(eps) == (32, 32)
            @test minimum(eps) ≈ 1.0
            @test maximum(eps) ≈ 9.0
        end

        @testset "1D Photonic" begin
            lat = lattice_1d(1.0)
            geo = Geometry(lat, Dielectric(1.0), [
                (Segment(0.0, 0.3), Dielectric(4.0))
            ])
            solver = Solver(Photonic1D(), geo, (64,); cutoff=5)

            eps = get_epsilon_field(solver)

            @test eps isa Vector{Float64}
            @test length(eps) == 64
            @test minimum(eps) ≈ 1.0
            @test maximum(eps) ≈ 4.0
        end
    end

    @testset "get_material_field" begin
        @testset "Phononic - density" begin
            lat = square_lattice(1.0)
            geo = Geometry(lat, IsotropicElastic(ρ=1000.0, λ=1e9, μ=1e9), [
                (Circle([0.5, 0.5], 0.2), IsotropicElastic(ρ=2000.0, λ=2e9, μ=2e9))
            ])
            solver = Solver(SHWave(), geo, (32, 32); cutoff=5)

            rho = get_material_field(solver, :ρ)

            @test rho isa Matrix{Float64}
            @test size(rho) == (32, 32)
            @test minimum(rho) ≈ 1000.0
            @test maximum(rho) ≈ 2000.0
        end
    end

    # =========================================================================
    # fix_phase
    # =========================================================================
    @testset "fix_phase" begin
        @testset "method=:max makes maximum point real" begin
            # Create complex field
            field = [1.0 + 2.0im, 3.0 + 1.0im, 0.5 - 0.5im]

            fixed = fix_phase(field; method=:max)

            # Maximum amplitude point (2nd) should be real
            max_idx = argmax(abs.(field))
            @test abs(imag(fixed[max_idx])) < 1e-10
            @test real(fixed[max_idx]) > 0  # positive real
        end

        @testset "amplitude is preserved" begin
            field = randn(ComplexF64, 10)

            fixed = fix_phase(field; method=:max)

            @test abs.(field) ≈ abs.(fixed)
        end

        @testset "method=:center" begin
            field = reshape(randn(ComplexF64, 9), 3, 3)

            fixed = fix_phase(field; method=:center)

            # Center point should be real
            @test abs(imag(fixed[2, 2])) < 1e-10
        end
    end

    # =========================================================================
    # field_energy
    # =========================================================================
    @testset "field_energy" begin
        @testset "2D TM energy is positive" begin
            lat = square_lattice(1.0)
            geo = Geometry(lat, Dielectric(1.0), [
                (Circle([0.5, 0.5], 0.2), Dielectric(9.0))
            ])
            solver = Solver(TMWave(), geo, (32, 32); cutoff=5)
            freqs, vecs = solve_at_k_with_vectors(solver, [0.0, 0.0], DenseMethod())

            energy = field_energy(solver, vecs[:, 1])

            @test energy isa Matrix{Float64}
            @test size(energy) == (32, 32)
            @test all(energy .>= 0)
            @test sum(energy) > 0
        end
    end

end  # @testset "Field Reconstruction"
