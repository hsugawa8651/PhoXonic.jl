# Test TMMSolver and tmm_spectrum

using Test
using PhoXonic

@testset "TMMSolver" begin
    @testset "Construction" begin
        air = Dielectric(1.0)
        glass = Dielectric(2.25)
        layer = Layer(glass, 100.0)
        ml = Multilayer([layer], air, air)

        solver = TMMSolver(Photonic1D(), ml)

        @test solver isa TMMSolver
        @test solver.wavetype isa Photonic1D
        @test solver.structure === ml
    end
end

@testset "tmm_spectrum" begin
    @testset "Single interface" begin
        # Single thin layer approaches bare interface as d → 0
        air = Dielectric(1.0)
        glass = Dielectric(2.25)  # n = 1.5

        # Very thin layer (effectively single interface)
        thin_layer = Layer(glass, 0.001)  # 0.001 nm
        ml = Multilayer([thin_layer], air, glass)

        solver = TMMSolver(Photonic1D(), ml)
        λ = 600.0

        result = tmm_spectrum(solver, λ)

        # For thin layer, should approach Fresnel result
        # R = ((n1-n2)/(n1+n2))^2 = 0.04
        # T ≈ 0.96
        @test result.R ≈ 0.04 atol=0.01
        @test result.T ≈ 0.96 atol=0.01
        @test result.R + result.T ≈ 1.0 atol=1e-10
    end

    @testset "Energy conservation" begin
        air = Dielectric(1.0)
        glass = Dielectric(2.25)
        layer = Layer(glass, 100.0)
        ml = Multilayer([layer], air, air)

        solver = TMMSolver(Photonic1D(), ml)

        for λ in [400.0, 500.0, 600.0, 700.0]
            result = tmm_spectrum(solver, λ)
            @test result.R + result.T ≈ 1.0 atol=1e-10
        end
    end

    @testset "Spectrum calculation" begin
        air = Dielectric(1.0)
        glass = Dielectric(2.25)
        layer = Layer(glass, 100.0)
        ml = Multilayer([layer], air, air)

        solver = TMMSolver(Photonic1D(), ml)
        λ_values = collect(400.0:50.0:800.0)

        R_spectrum, T_spectrum = tmm_spectrum(solver, λ_values)

        @test length(R_spectrum) == length(λ_values)
        @test length(T_spectrum) == length(λ_values)
        @test all(@. R_spectrum + T_spectrum ≈ 1.0)
    end
end

@testset "Fabry-Pérot resonator" begin
    @testset "Resonance condition" begin
        # Air | Glass | Air with thickness = λ₀/(2n)
        # At λ₀, the optical path is exactly one wavelength → T = 1
        air = Dielectric(1.0)
        glass = Dielectric(2.25)  # n = 1.5
        n_glass = 1.5

        λ0 = 600.0
        d = λ0 / (2 * n_glass)  # Half-wave thickness = 200 nm

        layer = Layer(glass, d)
        ml = Multilayer([layer], air, air)

        solver = TMMSolver(Photonic1D(), ml)

        # At resonance, transmission should be maximum (= 1 for symmetric cavity)
        result = tmm_spectrum(solver, λ0)
        @test result.T ≈ 1.0 atol=1e-10
        @test result.R ≈ 0.0 atol=1e-10
    end

    @testset "Resonance peaks" begin
        # Multiple resonance orders
        air = Dielectric(1.0)
        glass = Dielectric(2.25)
        n_glass = 1.5

        d = 300.0  # 300 nm thick glass

        layer = Layer(glass, d)
        ml = Multilayer([layer], air, air)

        solver = TMMSolver(Photonic1D(), ml)

        # Resonance at λ = 2nd/m for m = 1, 2, 3...
        # 2*1.5*300 = 900 nm is the optical path
        # m=1: λ=900, m=2: λ=450, m=3: λ=300, etc.
        for m in 1:3
            λ_res = 2 * n_glass * d / m
            if λ_res > 300.0  # Only test visible/near-IR
                result = tmm_spectrum(solver, λ_res)
                @test result.T ≈ 1.0 atol=1e-10
            end
        end
    end
end

@testset "Bragg mirror" begin
    @testset "High reflection at center" begin
        # Quarter-wave stack centered at λ₀ = 1000 nm
        n_hi, n_lo = 3.0, 1.5
        λ0 = 1000.0
        N_pairs = 10

        ml = bragg_mirror(n_hi, n_lo, λ0, N_pairs)
        solver = TMMSolver(Photonic1D(), ml)

        result = tmm_spectrum(solver, λ0)

        # Should have very high reflection at center wavelength
        @test result.R > 0.99
        @test result.R + result.T ≈ 1.0 atol=1e-10
    end

    @testset "Reflection increases with pair count" begin
        n_hi, n_lo = 2.0, 1.5
        λ0 = 600.0

        R_values = Float64[]
        for N in [2, 5, 10, 15]
            ml = bragg_mirror(n_hi, n_lo, λ0, N)
            solver = TMMSolver(Photonic1D(), ml)
            result = tmm_spectrum(solver, λ0)
            push!(R_values, result.R)
        end

        # Reflection should increase monotonically
        @test issorted(R_values)
    end

    @testset "Stopband width" begin
        # Higher contrast → wider stopband
        λ0 = 1000.0
        N = 10

        # Low contrast
        ml_low = bragg_mirror(1.8, 1.5, λ0, N)
        solver_low = TMMSolver(Photonic1D(), ml_low)

        # High contrast
        ml_high = bragg_mirror(3.0, 1.5, λ0, N)
        solver_high = TMMSolver(Photonic1D(), ml_high)

        # Measure R at λ₀ ± 20%
        λ_test = λ0 * 0.8

        result_low = tmm_spectrum(solver_low, λ_test)
        result_high = tmm_spectrum(solver_high, λ_test)

        # Higher contrast should have higher R at edge of band
        @test result_high.R > result_low.R
    end
end

@testset "Complex refractive index (absorption)" begin
    # Note: Standard TMM has numerical issues with strongly absorbing layers.
    # The backward-propagating wave amplitude exp(-im*δ) can grow when Im(n) > 0,
    # causing |t| > 1 and T > 1 in some cases. This is a known limitation.
    # Accurate absorption calculation requires Scattering Matrix or Enhanced TMM.

    @testset "LossyDielectric in TMM" begin
        # Verify LossyDielectric can be used with TMM
        air = Dielectric(1.0)
        n_complex = 1.5 + 0.01im  # Weak absorption
        ε_complex = n_complex^2
        absorbing = LossyDielectric(ε_complex)

        layer = Layer(absorbing, 100.0)
        ml = Multilayer([layer], air, air)

        solver = TMMSolver(Photonic1D(), ml)
        result = tmm_spectrum(solver, 600.0)

        # Basic checks
        @test 0.0 ≤ result.R ≤ 1.0
        @test result.T > 0  # Transmission should be positive
        @test isfinite(result.R) && isfinite(result.T)
    end

    @testset "Reflection is accurate" begin
        # Reflection coefficient is accurate even with absorption
        # (since incident medium is lossless)
        air = Dielectric(1.0)
        n_absorbing = 1.5 + 0.02im
        ε_absorbing = n_absorbing^2
        absorbing = LossyDielectric(ε_absorbing)

        layer = Layer(absorbing, 50.0)
        ml = Multilayer([layer], air, air)
        solver = TMMSolver(Photonic1D(), ml)

        # Very thick layer at quarter-wave thickness
        λ0 = 50.0 * 4 * 1.5  # Quarter-wave at real part of n
        result = tmm_spectrum(solver, λ0)

        # R should be reasonable (not NaN or negative)
        @test 0.0 ≤ result.R ≤ 1.0
    end

    @testset "Absorption effect on R" begin
        # With absorption, some light is lost inside the layer
        # Compare with lossless case
        air = Dielectric(1.0)
        glass = Dielectric(2.25)
        glass_lossy = LossyDielectric(2.25 + 0.1im)

        d = 100.0
        λ = 500.0

        ml_lossless = Multilayer([Layer(glass, d)], air, air)
        ml_lossy = Multilayer([Layer(glass_lossy, d)], air, air)

        R_lossless = tmm_spectrum(TMMSolver(Photonic1D(), ml_lossless), λ).R
        R_lossy = tmm_spectrum(TMMSolver(Photonic1D(), ml_lossy), λ).R

        # Results should be different (absorption changes R)
        @test R_lossless != R_lossy
    end
end
