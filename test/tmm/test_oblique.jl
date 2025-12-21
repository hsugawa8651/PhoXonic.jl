# TMM oblique incidence tests

using Test
using PhoXonic
using LinearAlgebra

@testset "Snell's law" begin
    @testset "Basic refraction" begin
        n1, n2 = 1.0, 1.5  # air to glass

        # θ1 = 30°
        θ1 = deg2rad(30)
        θ2 = PhoXonic.snell_angle(n1, n2, θ1)

        # Verify Snell's law: n1 sin θ1 = n2 sin θ2
        @test isapprox(n1 * sin(θ1), n2 * sin(real(θ2)), atol=1e-10)

        # θ2 should be smaller when entering denser medium (real part)
        @test real(θ2) < θ1

        # For propagating wave, imaginary part should be ~0
        @test abs(imag(θ2)) < 1e-10
    end

    @testset "Normal incidence" begin
        n1, n2 = 1.0, 1.5
        θ1 = 0.0
        θ2 = PhoXonic.snell_angle(n1, n2, θ1)

        @test real(θ2) ≈ 0.0 atol=1e-10
        @test abs(imag(θ2)) < 1e-10
    end

    @testset "Symmetric refraction" begin
        n1, n2 = 1.0, 1.5
        θ1 = deg2rad(30)

        # Air → Glass
        θ2 = PhoXonic.snell_angle(n1, n2, θ1)

        # Glass → Air (reverse)
        θ1_back = PhoXonic.snell_angle(n2, n1, θ2)

        @test isapprox(real(θ1_back), θ1, atol=1e-10)
    end
end

@testset "Total internal reflection" begin
    @testset "Critical angle" begin
        n1, n2 = 1.5, 1.0  # glass to air

        # Critical angle: θc = asin(n2/n1)
        θc = asin(n2 / n1)
        @test θc ≈ deg2rad(41.81) atol=0.01

        # At critical angle, θ2 = 90°
        θ2 = PhoXonic.snell_angle(n1, n2, θc)
        @test isapprox(real(θ2), π/2, atol=1e-6)
    end

    @testset "Above critical angle" begin
        n1, n2 = 1.5, 1.0
        θc = asin(n2 / n1)

        # Above critical angle: total internal reflection
        θ1 = θc + deg2rad(5)
        θ2 = PhoXonic.snell_angle(n1, n2, θ1)

        # θ2 should be complex (evanescent wave)
        @test imag(θ2) != 0.0
    end

    @testset "TIR reflectance" begin
        n1, n2 = 1.5, 1.0
        θc = asin(n2 / n1)
        θ1 = θc + deg2rad(10)  # above critical angle

        # Both TE and TM should have R = 1.0
        r_te, _ = PhoXonic.fresnel_oblique(n1, n2, θ1, :TE)
        r_tm, _ = PhoXonic.fresnel_oblique(n1, n2, θ1, :TM)

        @test abs(r_te)^2 ≈ 1.0 atol=1e-10
        @test abs(r_tm)^2 ≈ 1.0 atol=1e-10
    end
end

@testset "TE polarization (s-wave)" begin
    @testset "Normal incidence matches standard Fresnel" begin
        n1, n2 = 1.0, 1.5
        θ1 = 0.0

        r_te, t_te = PhoXonic.fresnel_oblique(n1, n2, θ1, :TE)
        r_normal, t_normal = PhoXonic.fresnel_coefficients(n1, n2)

        @test isapprox(r_te, r_normal, atol=1e-10)
        @test isapprox(t_te, t_normal, atol=1e-10)
    end

    @testset "Oblique incidence formulas" begin
        n1, n2 = 1.0, 1.5
        θ1 = deg2rad(45)
        θ2 = PhoXonic.snell_angle(n1, n2, θ1)

        r_te, t_te = PhoXonic.fresnel_oblique(n1, n2, θ1, :TE)

        # Theoretical: r_s = (n1 cos θ1 - n2 cos θ2) / (n1 cos θ1 + n2 cos θ2)
        r_theory = (n1 * cos(θ1) - n2 * cos(θ2)) / (n1 * cos(θ1) + n2 * cos(θ2))
        t_theory = 2 * n1 * cos(θ1) / (n1 * cos(θ1) + n2 * cos(θ2))

        @test isapprox(r_te, r_theory, atol=1e-10)
        @test isapprox(t_te, t_theory, atol=1e-10)
    end

    @testset "Energy conservation" begin
        n1, n2 = 1.0, 1.5

        for θ1_deg in [0, 15, 30, 45, 60, 75]
            θ1 = deg2rad(θ1_deg)
            θ2 = PhoXonic.snell_angle(n1, n2, θ1)

            r_te, t_te = PhoXonic.fresnel_oblique(n1, n2, θ1, :TE)

            R = abs(r_te)^2
            T = (n2 * cos(real(θ2))) / (n1 * cos(θ1)) * abs(t_te)^2

            @test isapprox(R + T, 1.0, atol=1e-10)
        end
    end
end

@testset "TM polarization (p-wave)" begin
    @testset "Normal incidence matches standard Fresnel" begin
        n1, n2 = 1.0, 1.5
        θ1 = 0.0

        r_tm, t_tm = PhoXonic.fresnel_oblique(n1, n2, θ1, :TM)
        r_normal, t_normal = PhoXonic.fresnel_coefficients(n1, n2)

        # TM and normal incidence have different sign conventions
        # but |r| and |t| should match, giving same R and T
        @test isapprox(abs(r_tm), abs(r_normal), atol=1e-10)
        @test isapprox(abs(t_tm), abs(t_normal), atol=1e-10)
    end

    @testset "Oblique incidence formulas" begin
        n1, n2 = 1.0, 1.5
        θ1 = deg2rad(45)
        θ2 = PhoXonic.snell_angle(n1, n2, θ1)

        r_tm, t_tm = PhoXonic.fresnel_oblique(n1, n2, θ1, :TM)

        # Theoretical: r_p = (n2 cos θ1 - n1 cos θ2) / (n2 cos θ1 + n1 cos θ2)
        r_theory = (n2 * cos(θ1) - n1 * cos(θ2)) / (n2 * cos(θ1) + n1 * cos(θ2))
        t_theory = 2 * n1 * cos(θ1) / (n2 * cos(θ1) + n1 * cos(θ2))

        @test isapprox(r_tm, r_theory, atol=1e-10)
        @test isapprox(t_tm, t_theory, atol=1e-10)
    end

    @testset "Brewster angle" begin
        n1, n2 = 1.0, 1.5

        # Brewster angle: θ_B = atan(n2/n1)
        θ_B = atan(n2 / n1)
        @test θ_B ≈ deg2rad(56.31) atol=0.01

        r_tm, _ = PhoXonic.fresnel_oblique(n1, n2, θ_B, :TM)

        # At Brewster angle, r_p = 0
        @test abs(r_tm) < 1e-10
    end

    @testset "TE vs TM at Brewster angle" begin
        n1, n2 = 1.0, 1.5
        θ_B = atan(n2 / n1)

        r_te, _ = PhoXonic.fresnel_oblique(n1, n2, θ_B, :TE)
        r_tm, _ = PhoXonic.fresnel_oblique(n1, n2, θ_B, :TM)

        # TM is zero, but TE is not
        @test abs(r_tm) < 1e-10
        @test abs(r_te) > 0.1
    end

    @testset "Energy conservation" begin
        n1, n2 = 1.0, 1.5

        for θ1_deg in [0, 15, 30, 45, 60, 75]
            θ1 = deg2rad(θ1_deg)
            θ2 = PhoXonic.snell_angle(n1, n2, θ1)

            r_tm, t_tm = PhoXonic.fresnel_oblique(n1, n2, θ1, :TM)

            R = abs(r_tm)^2
            T = (n2 * cos(real(θ2))) / (n1 * cos(θ1)) * abs(t_tm)^2

            @test isapprox(R + T, 1.0, atol=1e-10)
        end
    end
end

@testset "Oblique TMM solver" begin
    @testset "Single interface at oblique incidence" begin
        # Air → Glass interface
        n_air = 1.0
        n_glass = 1.5

        mat_air = Dielectric(n_air^2)
        mat_glass = Dielectric(n_glass^2)

        # Very thin layer (effectively single interface)
        ml = Multilayer([Layer(mat_glass, 1e-6)], mat_air, mat_glass)
        solver = TMMSolver(Photonic1D(), ml)

        θ1 = deg2rad(45)
        λ = 1.0

        # TE polarization
        result_te = tmm_spectrum(solver, λ; angle=θ1, polarization=:TE)
        r_te, _ = PhoXonic.fresnel_oblique(n_air, n_glass, θ1, :TE)
        @test isapprox(result_te.R, abs(r_te)^2, atol=0.01)

        # TM polarization
        result_tm = tmm_spectrum(solver, λ; angle=θ1, polarization=:TM)
        r_tm, _ = PhoXonic.fresnel_oblique(n_air, n_glass, θ1, :TM)
        @test isapprox(result_tm.R, abs(r_tm)^2, atol=0.01)
    end

    @testset "Brewster angle in TMM" begin
        n1, n2 = 1.0, 1.5
        θ_B = atan(n2 / n1)

        mat1 = Dielectric(n1^2)
        mat2 = Dielectric(n2^2)

        ml = Multilayer([Layer(mat2, 0.5)], mat1, mat2)
        solver = TMMSolver(Photonic1D(), ml)

        # At Brewster angle, TM reflectance should be ~0
        result_tm = tmm_spectrum(solver, 1.0; angle=θ_B, polarization=:TM)
        @test result_tm.R < 0.01

        # TE should still reflect
        result_te = tmm_spectrum(solver, 1.0; angle=θ_B, polarization=:TE)
        @test result_te.R > 0.01
    end

    @testset "Bragg mirror at oblique incidence" begin
        n_hi, n_lo = 2.5, 1.5
        λ0 = 1.0  # design wavelength at normal incidence

        mat_hi = Dielectric(n_hi^2)
        mat_lo = Dielectric(n_lo^2)

        # Quarter-wave thicknesses at normal incidence
        d_hi = λ0 / (4 * n_hi)
        d_lo = λ0 / (4 * n_lo)

        unit_cell = [Layer(mat_hi, d_hi), Layer(mat_lo, d_lo)]
        ml = periodic_multilayer(unit_cell, 10)
        solver = TMMSolver(Photonic1D(), ml)

        # At normal incidence, high reflection at λ0
        result_normal = tmm_spectrum(solver, λ0)
        @test result_normal.R > 0.99

        # At oblique incidence, stopband shifts to shorter wavelength
        θ = deg2rad(30)
        result_te = tmm_spectrum(solver, λ0; angle=θ, polarization=:TE)
        result_tm = tmm_spectrum(solver, λ0; angle=θ, polarization=:TM)

        # Reflectance should decrease at λ0 due to band shift
        @test result_te.R < result_normal.R || result_tm.R < result_normal.R
    end
end
