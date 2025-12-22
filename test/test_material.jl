# Material type tests

using Test
using PhoXonic

@testset "LossyDielectric" begin
    @testset "Construction" begin
        # Complex permittivity
        m1 = LossyDielectric(2.25 + 0.01im)
        @test m1 isa LossyDielectric
        @test m1.ε ≈ 2.25 + 0.01im
        @test m1.μ ≈ 1.0 + 0.0im

        # Real value converts to ComplexF64
        m2 = LossyDielectric(2.25)
        @test m2.ε isa ComplexF64
        @test m2.ε ≈ 2.25 + 0.0im

        # With permeability
        m3 = LossyDielectric(2.25 + 0.01im, 1.5)
        @test m3.μ ≈ 1.5 + 0.0im
    end

    @testset "Accessors" begin
        n = 1.5 + 0.1im
        ε = n^2  # ≈ 2.24 + 0.3im
        m = LossyDielectric(ε)

        @test permittivity(m) ≈ ε
        @test permeability(m) ≈ 1.0 + 0.0im
        @test refractive_index(m) ≈ n
    end

    @testset "Type conversion" begin
        # Convert Dielectric → LossyDielectric
        d = Dielectric(2.25)
        ld = convert(LossyDielectric, d)

        @test ld isa LossyDielectric
        @test ld.ε ≈ 2.25 + 0.0im
        @test ld.μ ≈ 1.0 + 0.0im
    end

    @testset "Type promotion" begin
        # promote_type(Dielectric, LossyDielectric) → LossyDielectric
        @test promote_type(Dielectric, LossyDielectric) == LossyDielectric
        @test promote_type(LossyDielectric, Dielectric) == LossyDielectric
    end
end

@testset "LossyDielectric with TMM Layer" begin
    @testset "Layer construction" begin
        m = LossyDielectric(2.25 + 0.01im)
        layer = Layer(m, 100.0)

        @test layer isa Layer{LossyDielectric}
        @test thickness(layer) == 100.0
        @test material(layer) isa LossyDielectric
    end

    @testset "Multilayer with mixed materials" begin
        # Mix Dielectric and LossyDielectric → should auto-promote
        air = Dielectric(1.0)
        glass = Dielectric(2.25)
        metal = LossyDielectric(-10.0 + 1.0im)

        layers = [Layer(glass, 100.0), Layer(metal, 20.0), Layer(glass, 100.0)]
        ml = Multilayer(layers, air, air)

        # All materials should be promoted to LossyDielectric
        @test ml isa Multilayer{LossyDielectric}
        @test ml.incident isa LossyDielectric
        @test ml.substrate isa LossyDielectric
    end
end

@testset "LossyDielectric restriction" begin
    @testset "PWE Solver should reject LossyDielectric" begin
        lat = square_lattice(1.0)
        lossy = LossyDielectric(2.25 + 0.1im)

        # Creating Geometry with LossyDielectric should fail for PWE
        @test_throws Union{ArgumentError,MethodError} begin
            geo = Geometry(lat, lossy)
            Solver(TEWave(), geo, (32, 32); cutoff=3)
        end
    end
end

@testset "Dielectric error fallbacks" begin
    @test_throws ErrorException Dielectric()  # no arguments
    @test_throws ErrorException Dielectric("invalid")
    @test_throws ErrorException Dielectric([1.0, 2.0])  # wrong type
end
