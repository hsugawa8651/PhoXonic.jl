# test/dimension_guards_test.jl
# The (dimension, wave type) lattice is not a free product: each wave type
# belongs to exactly one dimension, and each analysis function is defined only
# for the dimensions where the quantity exists. These tests pin that down.
#
# Run independently: julia --project=. test/dimension_guards_test.jl

using Test
using PhoXonic

@testset "dimension guards" begin
    steel = IsotropicElastic(; ρ=7800.0, λ=115e9, μ=82e9)

    lat1 = lattice_1d(1.0)
    geo1 = Geometry(lat1, Dielectric(1.0), [(Segment(0.2, 0.8), Dielectric(9.0))])
    geo1_el = Geometry(lat1, steel, [(Segment(0.2, 0.8), steel)])

    lat2 = square_lattice(1.0)
    geo2 = Geometry(lat2, Dielectric(1.0), [(Circle([0.5, 0.5], 0.2), Dielectric(9.0))])
    geo2_el = Geometry(lat2, steel, [(Circle([0.5, 0.5], 0.2), steel)])

    lat3 = cubic_lattice(1.0)
    geo3 = Geometry(
        lat3, Dielectric(1.0), [(Sphere([0.5, 0.5, 0.5], 0.2), Dielectric(9.0))]
    )
    geo3_el = Geometry(lat3, steel, [(Sphere([0.5, 0.5, 0.5], 0.2), steel)])

    @testset "the nine valid pairings construct" begin
        @test Solver(Photonic1D(), geo1, (64,); cutoff=5) isa Solver{Dim1,Photonic1D}
        @test Solver(Longitudinal1D(), geo1_el, (64,); cutoff=5) isa
            Solver{Dim1,Longitudinal1D}
        @test Solver(TEWave(), geo2, (16, 16); cutoff=4) isa Solver{Dim2,TEWave}
        @test Solver(TMWave(), geo2, (16, 16); cutoff=4) isa Solver{Dim2,TMWave}
        @test Solver(SHWave(), geo2_el, (16, 16); cutoff=4) isa Solver{Dim2,SHWave}
        @test Solver(PSVWave(), geo2_el, (16, 16); cutoff=4) isa Solver{Dim2,PSVWave}
        @test Solver(TransverseEM(), geo3, (8, 8, 8); cutoff=2) isa
            Solver{Dim3,TransverseEM}
        @test Solver(FullVectorEM(), geo3, (8, 8, 8); cutoff=2) isa
            Solver{Dim3,FullVectorEM}
        @test Solver(FullElastic(), geo3_el, (8, 8, 8); cutoff=2) isa
            Solver{Dim3,FullElastic}
    end

    @testset "a wave type outside its dimension is rejected at construction" begin
        # prepare_materials has no method for these pairs; its generic fallback
        # names both halves of the pairing.
        @test_throws ArgumentError Solver(TEWave(), geo1, (64,); cutoff=5)
        @test_throws ArgumentError Solver(SHWave(), geo1_el, (64,); cutoff=5)
        @test_throws ArgumentError Solver(Photonic1D(), geo2, (16, 16); cutoff=4)
        @test_throws ArgumentError Solver(TransverseEM(), geo2, (16, 16); cutoff=4)
        @test_throws ArgumentError Solver(TMWave(), geo3, (8, 8, 8); cutoff=2)
        @test_throws ArgumentError Solver(FullElastic(), geo2_el, (16, 16); cutoff=4)

        @test_throws "TEWave with Dim1" Solver(TEWave(), geo1, (64,); cutoff=5)
    end

    @testset "topological invariants are guarded by dispatch" begin
        s1 = Solver(Photonic1D(), geo1, (64,); cutoff=5)
        s2 = Solver(TMWave(), geo2, (16, 16); cutoff=4)
        s3 = Solver(TransverseEM(), geo3, (8, 8, 8); cutoff=2)

        # Zak phase exists only in 1D, Wilson loop only in 2D. A wrong dimension
        # is a MethodError at the call site, never a silently wrong number.
        @test_throws MethodError compute_zak_phase(s2, 1:2)
        @test_throws MethodError compute_zak_phase(s3, 1:2)
        @test_throws MethodError compute_wilson_spectrum(s1, 1:2)
        @test_throws MethodError compute_wilson_spectrum(s3, 1:2)
    end

    @testset "the Green's function is dimension-generic, the DOS is not" begin
        s3 = Solver(TransverseEM(), geo3, (8, 8, 8), DenseMethod(); cutoff=2)
        ω_values = range(0.5, 1.5; length=3)

        source = zeros(ComplexF64, matrix_dimension(s3))
        source[1] = 1.0
        G = compute_greens_function(s3, [0.1, 0.1, 0.1], ω_values, source)
        @test length(G) == length(ω_values)
        @test all(isfinite, abs.(G[1]))

        # The density of states has no 3D method, in core or in the extension.
        k_points = [[0.1, 0.1, 0.1]]
        @test_throws MethodError compute_dos(s3, ω_values, k_points)
        @test_throws MethodError compute_dos_stochastic(s3, ω_values, k_points)
        @test_throws MethodError compute_ldos(s3, [0.5, 0.5, 0.5], ω_values, k_points)
    end
end
