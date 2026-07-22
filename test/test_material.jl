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

@testset "MultiphysicsMaterial" begin
    # IsotropicElastic positional order is (rho, C11, C12, C44)
    si_ph = Dielectric(11.7)
    si_el = IsotropicElastic(2330.0, 165.7e9, 63.9e9, 79.6e9)

    @testset "Construction" begin
        # Both argument orders construct, and both store the same way
        mp = MultiphysicsMaterial(si_ph, si_el)
        pm = MultiphysicsMaterial(si_el, si_ph)

        @test mp isa MultiphysicsMaterial
        @test pm isa MultiphysicsMaterial
        @test mp isa Material

        for m in (mp, pm)
            @test m.photonic isa PhotonicMaterial
            @test m.elastic isa ElasticMaterial
            @test m.photonic === si_ph
            @test m.elastic === si_el
        end

        # A void is a legitimate elastic side
        @test MultiphysicsMaterial(Dielectric(1.0), ElasticVoid()) isa MultiphysicsMaterial
    end

    @testset "Same-family pairs do not construct" begin
        @test_throws MethodError MultiphysicsMaterial(Dielectric(1.0), Dielectric(11.7))
        @test_throws MethodError MultiphysicsMaterial(si_el, ElasticVoid())
    end

    @testset "Property delegation" begin
        mp = MultiphysicsMaterial(si_ph, si_el)

        # Photonic side: the four symbols discretize actually asks for
        @test PhoXonic.get_property(mp, :ε) == PhoXonic.get_property(si_ph, :ε)
        @test PhoXonic.get_property(mp, :μ) == PhoXonic.get_property(si_ph, :μ)
        @test PhoXonic.get_property(mp, :ε_inv) == PhoXonic.get_property(si_ph, :ε_inv)
        @test PhoXonic.get_property(mp, :μ_inv) == PhoXonic.get_property(si_ph, :μ_inv)

        # Elastic side: the four symbols discretize actually asks for
        @test PhoXonic.get_property(mp, :ρ) == PhoXonic.get_property(si_el, :ρ)
        @test PhoXonic.get_property(mp, :C11) == PhoXonic.get_property(si_el, :C11)
        @test PhoXonic.get_property(mp, :C12) == PhoXonic.get_property(si_el, :C12)
        @test PhoXonic.get_property(mp, :C44) == PhoXonic.get_property(si_el, :C44)

        # An unknown property must still be an error, not a silent fallback
        @test_throws Union{ErrorException,ArgumentError} PhoXonic.get_property(mp, :nope)
    end
end

@testset "Geometry material-class validation" begin
    lat = square_lattice(1.0)
    circle = Circle([0.0, 0.0], 0.2)

    air = Dielectric(1.0)
    rod = Dielectric(11.7)
    lossy = LossyDielectric(2.25 + 0.01im)

    steel = IsotropicElastic(7800.0, 280e9, 130e9, 75e9)
    void = ElasticVoid()

    mp_bg = MultiphysicsMaterial(air, void)
    mp_incl = MultiphysicsMaterial(rod, steel)

    @testset "Coherent geometries are accepted" begin
        # All photonic
        @test Geometry(lat, air, [(circle, rod)]) isa Geometry
        # Photonic promotion across Dielectric and LossyDielectric stays legal
        @test Geometry(lat, air, [(circle, lossy)]) isa Geometry
        # All elastic, including the void
        @test Geometry(lat, steel, [(circle, steel)]) isa Geometry
        @test Geometry(lat, steel, [(circle, void)]) isa Geometry
        # All multiphysics
        @test Geometry(lat, mp_bg, [(circle, mp_incl)]) isa Geometry
    end

    @testset "Mixed classes are rejected at construction" begin
        @test_throws ArgumentError Geometry(lat, air, [(circle, steel)])
        @test_throws ArgumentError Geometry(lat, steel, [(circle, air)])
        @test_throws ArgumentError Geometry(lat, mp_bg, [(circle, rod)])
        @test_throws ArgumentError Geometry(lat, air, [(circle, mp_incl)])
        @test_throws ArgumentError Geometry(lat, mp_bg, [(circle, steel)])
        @test_throws ArgumentError Geometry(lat, steel, [(circle, mp_incl)])
    end

    @testset "The message names the problem" begin
        err = try
            Geometry(lat, air, [(circle, steel)])
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("Cannot mix material classes in one Geometry", err.msg)
    end
end

@testset "Solver wave/material compatibility" begin
    lat = square_lattice(1.0)
    circle = Circle([0.0, 0.0], 0.2)
    res = (16, 16)

    air = Dielectric(1.0)
    rod = Dielectric(11.7)
    steel = IsotropicElastic(7800.0, 280e9, 130e9, 75e9)
    void = ElasticVoid()

    geo_ph = Geometry(lat, air, [(circle, rod)])
    geo_el = Geometry(lat, steel, [(circle, void)])
    geo_mp = Geometry(
        lat, MultiphysicsMaterial(air, void), [(circle, MultiphysicsMaterial(rod, steel))]
    )

    @testset "Matching families construct" begin
        @test Solver(TEWave(), geo_ph, res; cutoff=3) isa Solver
        @test Solver(TMWave(), geo_ph, res; cutoff=3) isa Solver
        @test Solver(SHWave(), geo_el, res; cutoff=3) isa Solver
        @test Solver(PSVWave(), geo_el, res; cutoff=3) isa Solver
    end

    @testset "Multiphysics satisfies either family" begin
        @test Solver(TEWave(), geo_mp, res; cutoff=3) isa Solver
        @test Solver(TMWave(), geo_mp, res; cutoff=3) isa Solver
        @test Solver(SHWave(), geo_mp, res; cutoff=3) isa Solver
        @test Solver(PSVWave(), geo_mp, res; cutoff=3) isa Solver
    end

    @testset "Mismatched families are rejected at Solver construction" begin
        @test_throws ArgumentError Solver(SHWave(), geo_ph, res; cutoff=3)
        @test_throws ArgumentError Solver(PSVWave(), geo_ph, res; cutoff=3)
        @test_throws ArgumentError Solver(TEWave(), geo_el, res; cutoff=3)
        @test_throws ArgumentError Solver(TMWave(), geo_el, res; cutoff=3)
    end

    @testset "The message names the wave and the structure" begin
        err = try
            Solver(SHWave(), geo_ph, res; cutoff=3)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("SHWave", err.msg)
        @test occursin("elastic", err.msg)
    end
end
