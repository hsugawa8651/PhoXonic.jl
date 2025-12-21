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
