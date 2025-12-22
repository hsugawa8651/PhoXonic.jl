# Error fallback tests
# Run independently: julia --project=. test/test_error_fallbacks.jl

using Test
using PhoXonic

@testset "Error Fallbacks" begin

    @testset "Shape dimension mismatch" begin
        # Verify that mismatched dimensions throw ArgumentError with helpful message
        air = Dielectric(1.0)

        @testset "2D lattice rejects 3D shapes" begin
            lat2d = square_lattice(1.0)
            @test_throws ArgumentError Geometry(lat2d, air, [(Sphere([0,0,0], 0.2), air)])
            @test_throws ArgumentError Geometry(lat2d, air, [(Cylinder([0,0,0], 0.2, 0.5), air)])
            @test_throws ArgumentError Geometry(lat2d, air, [(Slab(0.0, 0.5), air)])
        end

        @testset "3D lattice rejects 2D shapes" begin
            lat3d = cubic_lattice(1.0)
            @test_throws ArgumentError Geometry(lat3d, air, [(Circle([0,0], 0.2), air)])
            @test_throws ArgumentError Geometry(lat3d, air, [(Rectangle([0,0], [0.3,0.3]), air)])
        end

        @testset "1D lattice rejects 2D/3D shapes" begin
            lat1d = lattice_1d(1.0)
            @test_throws ArgumentError Geometry(lat1d, air, [(Circle([0,0], 0.2), air)])
            @test_throws ArgumentError Geometry(lat1d, air, [(Sphere([0,0,0], 0.2), air)])
        end

        @testset "2D/3D lattice rejects 1D shapes" begin
            lat2d = square_lattice(1.0)
            lat3d = cubic_lattice(1.0)
            @test_throws ArgumentError Geometry(lat2d, air, [(Segment(0.0, 0.5), air)])
            @test_throws ArgumentError Geometry(lat3d, air, [(Segment(0.0, 0.5), air)])
        end
    end

    @testset "GFMethod fallback" begin
        # Define a custom unsupported method for testing
        struct UnsupportedGF <: GFMethod end

        # Setup minimal solver
        lat = square_lattice(1.0)
        air = Dielectric(1.0)
        geo = Geometry(lat, air)
        solver = Solver(TEWave(), geo, (8, 8))
        k = [0.0, 0.0]
        ω_values = [0.5, 1.0]
        source = ones(ComplexF64, solver.basis.num_pw)
        position = [0.5, 0.5]
        k_points = [[0.0, 0.0], [0.5, 0.0]]

        @test_throws ArgumentError compute_greens_function(solver, k, ω_values, source, UnsupportedGF())
        @test_throws ArgumentError compute_dos(solver, ω_values, k_points, UnsupportedGF())
        @test_throws ArgumentError compute_ldos(solver, position, ω_values, k_points, UnsupportedGF())
    end

    @testset "High-level API fallbacks" begin
        @testset "Solver constructor" begin
            @test_throws ErrorException Solver("not a wave", "not a geometry", (64, 64))
            @test_throws ErrorException Solver(123, nothing, (64, 64))
            @test_throws ErrorException Solver(:invalid, [], 64)
        end

        @testset "compute_bands" begin
            @test_throws ErrorException compute_bands("not a solver", "not a kpath")
            @test_throws ErrorException compute_bands(123, [0.0, 0.0])
        end

        @testset "discretize" begin
            @test_throws ErrorException discretize("not a geometry", "not a basis")
            @test_throws ErrorException discretize(nothing, nothing, nothing)
        end
    end

    @testset "Mid-level API fallbacks" begin
        @testset "solve_at_k" begin
            @test_throws ErrorException solve_at_k("not a solver", [0.0, 0.0], DenseMethod())
            @test_throws ErrorException solve_at_k(123, [0.0, 0.0], DenseMethod())
            @test_throws ErrorException solve_at_k(nothing, 0.0, nothing)
        end

        @testset "matrix_dimension" begin
            @test_throws ErrorException matrix_dimension("not a solver")
            @test_throws ErrorException matrix_dimension(123)
            @test_throws ErrorException matrix_dimension([1, 2, 3])
        end
    end

end
