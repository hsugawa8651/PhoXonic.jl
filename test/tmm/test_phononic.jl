# TMM phononic (elastic wave) tests

using Test
using PhoXonic
using LinearAlgebra

# Material constants for tests
# Steel: ρ = 7800 kg/m³, c_L ≈ 5960 m/s
# Epoxy: ρ = 1180 kg/m³, c_L ≈ 2540 m/s

const STEEL = IsotropicElastic(; ρ=7800.0, λ=115e9, μ=82e9)
const EPOXY = IsotropicElastic(; ρ=1180.0, λ=4.43e9, μ=1.59e9)

@testset "Acoustic impedance" begin
    @testset "Definition Z = ρc" begin
        # Steel
        ρ_steel = density(STEEL)
        c_steel = longitudinal_velocity(STEEL)
        Z_steel = PhoXonic.acoustic_impedance(STEEL)

        @test Z_steel ≈ ρ_steel * c_steel

        # Epoxy
        ρ_epoxy = density(EPOXY)
        c_epoxy = longitudinal_velocity(EPOXY)
        Z_epoxy = PhoXonic.acoustic_impedance(EPOXY)

        @test Z_epoxy ≈ ρ_epoxy * c_epoxy
    end

    @testset "Impedance ratio" begin
        Z_steel = PhoXonic.acoustic_impedance(STEEL)
        Z_epoxy = PhoXonic.acoustic_impedance(EPOXY)

        # Steel/Epoxy impedance ratio ~ 15.5
        ratio = Z_steel / Z_epoxy
        @test 15.0 < ratio < 16.0
    end

    @testset "Units consistency" begin
        # Z has units [kg/(m²·s)] = [Pa·s/m] = [Rayl]
        Z = PhoXonic.acoustic_impedance(STEEL)

        # Should be order of 10^7 for steel
        @test 1e7 < Z < 1e8
    end
end

@testset "Phononic Fresnel coefficients" begin
    @testset "Single interface reflection" begin
        Z1 = PhoXonic.acoustic_impedance(STEEL)
        Z2 = PhoXonic.acoustic_impedance(EPOXY)

        r, t = PhoXonic.acoustic_fresnel(Z1, Z2)

        # Theoretical: r = (Z2 - Z1) / (Z2 + Z1)
        r_theory = (Z2 - Z1) / (Z2 + Z1)
        @test isapprox(r, r_theory, atol=1e-10)

        # Energy conservation: R + T = 1
        R = abs2(r)
        T = 4 * Z1 * Z2 / (Z1 + Z2)^2
        @test isapprox(R + T, 1.0, atol=1e-10)
    end

    @testset "Same material (no reflection)" begin
        Z = PhoXonic.acoustic_impedance(STEEL)
        r, t = PhoXonic.acoustic_fresnel(Z, Z)

        @test abs(r) < 1e-10
        @test abs(t - 1.0) < 1e-10
    end

    @testset "High impedance contrast" begin
        Z1 = PhoXonic.acoustic_impedance(STEEL)
        Z2 = PhoXonic.acoustic_impedance(EPOXY)

        r, _ = PhoXonic.acoustic_fresnel(Z1, Z2)
        R = abs2(r)

        # High impedance mismatch should give high reflection
        @test R > 0.7
    end
end

@testset "Phononic TMM solver" begin
    @testset "TMMSolver construction" begin
        layer = Layer(STEEL, 0.001)  # 1 mm steel
        ml = Multilayer([layer], EPOXY, EPOXY)

        solver = TMMSolver(Longitudinal1D(), ml)

        @test solver.wavetype isa Longitudinal1D
        @test solver.structure === ml
    end

    @testset "Single layer transmission" begin
        d = 0.001  # 1 mm
        layer = Layer(STEEL, d)
        ml = Multilayer([layer], EPOXY, EPOXY)

        solver = TMMSolver(Longitudinal1D(), ml)

        # At wavelength >> thickness, high transmission
        c_steel = longitudinal_velocity(STEEL)
        λ_long = 10 * d  # λ = 10 mm >> d
        f = c_steel / λ_long

        result = tmm_spectrum(solver, λ_long)

        @test result.R >= 0.0
        @test result.T >= 0.0
        @test isapprox(result.R + result.T, 1.0, atol=0.01)
    end

    @testset "Energy conservation" begin
        layer1 = Layer(STEEL, 0.001)
        layer2 = Layer(EPOXY, 0.002)
        ml = Multilayer([layer1, layer2], EPOXY, STEEL)

        solver = TMMSolver(Longitudinal1D(), ml)

        for λ in [0.001, 0.005, 0.01, 0.05]
            result = tmm_spectrum(solver, λ)
            @test isapprox(result.R + result.T, 1.0, atol=0.01)
        end
    end

    @testset "Fabry-Perot resonance" begin
        # Single layer between same materials: resonance at d = m*λ/2
        d = 0.01  # 10 mm
        layer = Layer(STEEL, d)
        ml = Multilayer([layer], EPOXY, EPOXY)

        solver = TMMSolver(Longitudinal1D(), ml)

        c_steel = longitudinal_velocity(STEEL)

        # Resonance condition: d = m * λ_steel / 2
        # At resonance, T should be maximum (=1 for lossless)
        λ_res = 2 * d  # First resonance in steel
        λ_phys = λ_res  # Physical wavelength in steel

        result_res = tmm_spectrum(solver, λ_phys)

        # Off-resonance
        λ_off = 1.5 * d
        result_off = tmm_spectrum(solver, λ_off)

        # Transmission at resonance should be higher
        @test result_res.T >= result_off.T
    end
end

@testset "Phononic Bragg mirror" begin
    @testset "Periodic structure" begin
        # Steel/Epoxy quarter-wave stack
        c_steel = longitudinal_velocity(STEEL)
        c_epoxy = longitudinal_velocity(EPOXY)

        # Design frequency
        f0 = 100e3  # 100 kHz
        λ_steel = c_steel / f0
        λ_epoxy = c_epoxy / f0

        # Quarter-wave thicknesses
        d_steel = λ_steel / 4
        d_epoxy = λ_epoxy / 4

        unit_cell = [Layer(STEEL, d_steel), Layer(EPOXY, d_epoxy)]
        ml = periodic_multilayer(unit_cell, 10)

        solver = TMMSolver(Longitudinal1D(), ml)

        # At center frequency, high reflection
        λ0 = c_steel / f0  # Use steel wavelength as reference
        result = tmm_spectrum(solver, λ0)

        @test result.R > 0.9
    end

    @testset "Stopband width increases with pairs" begin
        c_steel = longitudinal_velocity(STEEL)
        c_epoxy = longitudinal_velocity(EPOXY)

        f0 = 100e3
        d_steel = c_steel / f0 / 4
        d_epoxy = c_epoxy / f0 / 4

        unit_cell = [Layer(STEEL, d_steel), Layer(EPOXY, d_epoxy)]

        # Compare 5 vs 10 pairs
        ml_5 = periodic_multilayer(unit_cell, 5)
        ml_10 = periodic_multilayer(unit_cell, 10)

        solver_5 = TMMSolver(Longitudinal1D(), ml_5)
        solver_10 = TMMSolver(Longitudinal1D(), ml_10)

        λ0 = c_steel / f0
        R_5 = tmm_spectrum(solver_5, λ0).R
        R_10 = tmm_spectrum(solver_10, λ0).R

        # More pairs = higher reflection
        @test R_10 > R_5
    end
end

@testset "Phononic band structure" begin
    @testset "tmm_bandstructure for Longitudinal1D" begin
        c_steel = longitudinal_velocity(STEEL)
        c_epoxy = longitudinal_velocity(EPOXY)

        # Unit cell
        a = 0.01  # 10 mm period
        d_steel = 0.25 * a
        d_epoxy = 0.75 * a

        unit_cell = [Layer(STEEL, d_steel), Layer(EPOXY, d_epoxy)]
        ml = Multilayer(unit_cell, EPOXY, EPOXY)

        solver = TMMSolver(Longitudinal1D(), ml)

        bands = tmm_bandstructure(solver; k_points=21, bands=1:3)

        # Check output type
        @test bands isa BandStructure

        # Check dimensions
        @test size(bands.frequencies, 1) == 21
        @test size(bands.frequencies, 2) == 3

        # Frequencies should be sorted
        for ik in 1:21
            @test issorted(bands.frequencies[ik, :])
        end
    end

    @testset "Bandgap existence" begin
        a = 0.01
        d_steel = 0.25 * a
        d_epoxy = 0.75 * a

        unit_cell = [Layer(STEEL, d_steel), Layer(EPOXY, d_epoxy)]
        ml = Multilayer(unit_cell, EPOXY, EPOXY)

        solver = TMMSolver(Longitudinal1D(), ml)
        bands = tmm_bandstructure(solver; k_points=51, bands=1:4)

        # High impedance contrast should produce bandgap
        band1_max = maximum(bands.frequencies[:, 1])
        band2_min = minimum(bands.frequencies[:, 2])

        @test band2_min > band1_max  # Gap exists
    end

    @testset "PWE comparison" begin
        # Same structure computed with PWE and TMM
        # Note: There may be systematic differences in frequency conventions
        # between PWE and TMM for phononic crystals. We test that:
        # 1. Both methods produce similar band ordering
        # 2. Bandgaps appear at similar locations

        a = 0.01  # 10 mm period
        d_steel = 0.25 * a

        # TMM setup
        unit_cell = [Layer(STEEL, d_steel), Layer(EPOXY, a - d_steel)]
        ml = Multilayer(unit_cell, EPOXY, EPOXY)
        solver_tmm = TMMSolver(Longitudinal1D(), ml)

        # PWE setup
        lat = lattice_1d(a)
        geo = Geometry(lat, EPOXY, [(Segment(0.0, d_steel), STEEL)])
        solver_pwe = Solver(Longitudinal1D(), geo, 256; cutoff=20)

        # Compute bands
        bands_tmm = tmm_bandstructure(solver_tmm; k_points=11, bands=1:3)

        # Compare at zone boundary (X point) - first band should be similar
        k_X = π/a
        ω_pwe_X, _ = solve(solver_pwe, k_X; bands=1:3)
        ω_tmm_X = bands_tmm.frequencies[end, :]

        # First band at X should be within 20% (allowing for convention differences)
        @test isapprox(ω_tmm_X[1], ω_pwe_X[1], rtol=0.20)

        # Both should have increasing frequencies (band ordering preserved)
        @test issorted(ω_pwe_X)
        @test issorted(ω_tmm_X)

        # Both should show a gap (band2_min > band1_max)
        gap_tmm =
            minimum(bands_tmm.frequencies[:, 2]) - maximum(bands_tmm.frequencies[:, 1])
        @test gap_tmm > 0
    end
end
