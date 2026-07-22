# Test Layer and Multilayer types

using Test
using PhoXonic

@testset "Layer" begin
    @testset "Photonic Layer" begin
        mat = Dielectric(9.0)
        layer = Layer(mat, 0.25)

        @test thickness(layer) == 0.25
        @test material(layer) === mat
        @test material(layer) isa Dielectric
    end

    @testset "Phononic Layer" begin
        mat = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)
        layer = Layer(mat, 0.5)

        @test thickness(layer) == 0.5
        @test material(layer) === mat
        @test material(layer) isa IsotropicElastic
    end
end

@testset "Multilayer" begin
    air = Dielectric(1.0)
    glass = Dielectric(2.25)

    layer1 = Layer(glass, 0.25)
    layer2 = Layer(air, 0.75)

    ml = Multilayer([layer1, layer2], air, air)

    @test nlayers(ml) == 2
    @test incident(ml) === air
    @test substrate(ml) === air
    @test total_thickness(ml) ≈ 1.0
    @test layers(ml)[1] === layer1
    @test layers(ml)[2] === layer2
end

@testset "Multilayer with mixed elastic materials" begin
    # A solid layer in a void: the ordinary phononic stack, and the mixture
    # promote_rule has always allowed. It lands on the abstract ElasticMaterial,
    # which is where the layer element type used to be lost.
    steel = IsotropicElastic(7800.0, 280e9, 130e9, 75e9)
    ml = Multilayer([Layer(steel, 1.0)], ElasticVoid(), ElasticVoid())

    @test ml isa Multilayer{ElasticMaterial}
    @test nlayers(ml) == 1
end

@testset "Helper functions" begin
    @testset "periodic_multilayer" begin
        air = Dielectric(1.0)
        glass = Dielectric(2.25)
        unit_cell = [Layer(glass, 0.25), Layer(air, 0.75)]

        ml = periodic_multilayer(unit_cell, 3)

        @test nlayers(ml) == 6  # 2 layers × 3 periods
        @test total_thickness(ml) ≈ 3.0
    end

    @testset "bragg_mirror" begin
        ml = bragg_mirror(3.0, 1.5, 1000.0, 5)

        # Quarter-wave condition: n*d = λ/4
        @test nlayers(ml) == 10  # 2 layers × 5 pairs

        # Check thicknesses: d = λ/(4n)
        d_hi = 1000.0 / (4 * 3.0)   # ≈ 83.33
        d_lo = 1000.0 / (4 * 1.5)   # ≈ 166.67

        @test thickness(layers(ml)[1]) ≈ d_hi
        @test thickness(layers(ml)[2]) ≈ d_lo
    end
end

@testset "Multilayer material-class validation" begin
    air = Dielectric(1.0)
    glass = Dielectric(2.25)
    steel = IsotropicElastic(7800.0, 280e9, 130e9, 75e9)
    epoxy = IsotropicElastic(1180.0, 7.61e9, 4.43e9, 1.59e9)

    mp_air = MultiphysicsMaterial(air, ElasticVoid())
    mp_glass = MultiphysicsMaterial(glass, steel)

    @testset "Coherent stacks are accepted" begin
        @test Multilayer([Layer(glass, 100.0)], air, air) isa Multilayer
        @test Multilayer([Layer(steel, 1.0)], epoxy, epoxy) isa Multilayer
        @test Multilayer([Layer(mp_glass, 100.0)], mp_air, mp_air) isa Multilayer
    end

    @testset "Mixed classes are rejected" begin
        @test_throws ArgumentError Multilayer([Layer(steel, 1.0)], air, air)
        @test_throws ArgumentError Multilayer([Layer(glass, 100.0)], epoxy, epoxy)
        @test_throws ArgumentError Multilayer([Layer(mp_glass, 100.0)], air, air)
    end
end

@testset "TMMSolver wave/material compatibility" begin
    air = Dielectric(1.0)
    glass = Dielectric(2.25)
    steel = IsotropicElastic(7800.0, 280e9, 130e9, 75e9)
    epoxy = IsotropicElastic(1180.0, 7.61e9, 4.43e9, 1.59e9)

    ml_ph = Multilayer([Layer(glass, 100.0)], air, air)
    ml_el = Multilayer([Layer(steel, 1.0)], epoxy, epoxy)
    ml_mp = Multilayer(
        [Layer(MultiphysicsMaterial(glass, steel), 100.0)],
        MultiphysicsMaterial(air, ElasticVoid()),
        MultiphysicsMaterial(air, ElasticVoid()),
    )

    @testset "Matching families construct" begin
        @test TMMSolver(Photonic1D(), ml_ph) isa TMMSolver
        @test TMMSolver(Longitudinal1D(), ml_el) isa TMMSolver
    end

    @testset "Multiphysics serves either family" begin
        @test TMMSolver(Photonic1D(), ml_mp) isa TMMSolver
        @test TMMSolver(Longitudinal1D(), ml_mp) isa TMMSolver
    end

    @testset "Mismatched families are rejected" begin
        @test_throws ArgumentError TMMSolver(Photonic1D(), ml_el)
        @test_throws ArgumentError TMMSolver(Longitudinal1D(), ml_ph)
    end

    @testset "Multiphysics actually computes, not just constructs" begin
        # Constructing is not the same as reaching the material. The photonic
        # side reads mat.ε and the elastic side reads the velocities, so both
        # paths have to be run, not merely built.
        solver = TMMSolver(Photonic1D(), ml_mp)
        r = tmm_spectrum(solver, 600.0)
        @test r.R >= 0.0
        @test r.T >= 0.0
        @test r.R + r.T + r.A ≈ 1.0 atol = 1e-8
    end
end

@testset "Multiphysics through the phononic TMM path" begin
    # Same materials and geometry as test_phononic.jl, wrapped so that each
    # region also carries an electromagnetic side. The elastic numbers must not
    # depend on that wrapping.
    steel = IsotropicElastic(; ρ=7800.0, λ=115e9, μ=82e9)
    epoxy = IsotropicElastic(; ρ=1180.0, λ=4.43e9, μ=1.59e9)

    d = 0.001          # 1 mm
    λ_long = 10 * d    # λ = 10 mm >> d

    mp_steel = MultiphysicsMaterial(Dielectric(11.7), steel)
    mp_epoxy = MultiphysicsMaterial(Dielectric(3.6), epoxy)

    ml_plain = Multilayer([Layer(steel, d)], epoxy, epoxy)
    ml_mp = Multilayer([Layer(mp_steel, d)], mp_epoxy, mp_epoxy)

    @testset "tmm_spectrum" begin
        plain = tmm_spectrum(TMMSolver(Longitudinal1D(), ml_plain), λ_long)
        mp = tmm_spectrum(TMMSolver(Longitudinal1D(), ml_mp), λ_long)

        @test mp.R ≈ plain.R
        @test mp.T ≈ plain.T
        @test mp.R + mp.T + mp.A ≈ 1.0 atol = 1e-8
    end

    @testset "tmm_bandstructure" begin
        # A separate entry point, and it reaches the material through its own
        # code path rather than through system_matrix_acoustic.
        unit = [Layer(steel, d), Layer(epoxy, 2d)]
        unit_mp = [Layer(mp_steel, d), Layer(mp_epoxy, 2d)]

        plain = tmm_bandstructure(
            TMMSolver(Longitudinal1D(), Multilayer(unit, epoxy, epoxy));
            k_points=5,
            bands=1:2,
        )
        mp = tmm_bandstructure(
            TMMSolver(Longitudinal1D(), Multilayer(unit_mp, mp_epoxy, mp_epoxy));
            k_points=5,
            bands=1:2,
        )

        @test mp.frequencies ≈ plain.frequencies
    end
end
