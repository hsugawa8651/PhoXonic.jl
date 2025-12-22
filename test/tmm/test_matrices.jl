# Test transfer matrix calculations

using Test
using PhoXonic
using LinearAlgebra

@testset "Fresnel coefficients" begin
    @testset "Air to Glass (n=1.5)" begin
        n1, n2 = 1.0, 1.5
        r, t = PhoXonic.fresnel_coefficients(n1, n2)

        # Theoretical: r = (n1 - n2)/(n1 + n2) = -0.2
        #              t = 2*n1/(n1 + n2) = 0.8
        @test r ≈ -0.2
        @test t ≈ 0.8

        # Power conservation: R + T = 1
        R = abs2(r)
        T = (n2/n1) * abs2(t)
        @test R + T ≈ 1.0
    end

    @testset "Glass to Air" begin
        n1, n2 = 1.5, 1.0
        r, t = PhoXonic.fresnel_coefficients(n1, n2)

        # r = (1.5 - 1)/(1.5 + 1) = 0.2
        # t = 2*1.5/(1.5 + 1) = 1.2
        @test r ≈ 0.2
        @test t ≈ 1.2

        # Power conservation
        R = abs2(r)
        T = (n2/n1) * abs2(t)
        @test R + T ≈ 1.0
    end

    @testset "High contrast" begin
        n1, n2 = 1.0, 3.0
        r, t = PhoXonic.fresnel_coefficients(n1, n2)

        # r = (1 - 3)/(1 + 3) = -0.5
        @test r ≈ -0.5

        # Power conservation
        R = abs2(r)
        T = (n2/n1) * abs2(t)
        @test R + T ≈ 1.0
    end
end

@testset "Propagation matrix" begin
    @testset "Phase delay" begin
        n = 1.5
        d = 100.0  # thickness
        λ = 600.0  # wavelength

        P = PhoXonic.propagation_matrix(n, d, λ)

        # Phase delay: δ = 2πnd/λ
        δ = 2π * n * d / λ

        # P = [exp(iδ) 0; 0 exp(-iδ)]
        @test P[1, 1] ≈ exp(im * δ)
        @test P[2, 2] ≈ exp(-im * δ)
        @test abs(P[1, 2]) < 1e-15
        @test abs(P[2, 1]) < 1e-15
    end

    @testset "Round trip identity" begin
        n = 2.0
        d = 50.0
        λ = 400.0

        P_forward = PhoXonic.propagation_matrix(n, d, λ)
        P_backward = PhoXonic.propagation_matrix(n, -d, λ)

        # P(d) * P(-d) = I
        @test P_forward * P_backward ≈ I
    end

    @testset "Zero thickness" begin
        n = 1.5
        d = 0.0
        λ = 500.0

        P = PhoXonic.propagation_matrix(n, d, λ)
        @test P ≈ I
    end
end

@testset "Interface matrix" begin
    @testset "Basic form" begin
        n1, n2 = 1.0, 1.5
        M = PhoXonic.interface_matrix(n1, n2)

        # M = (1/t) * [1 r; r 1] where r, t are Fresnel coeffs
        r, t = PhoXonic.fresnel_coefficients(n1, n2)

        @test M[1, 1] ≈ 1/t
        @test M[1, 2] ≈ r/t
        @test M[2, 1] ≈ r/t
        @test M[2, 2] ≈ 1/t
    end

    @testset "Determinant" begin
        n1, n2 = 1.0, 2.0
        M = PhoXonic.interface_matrix(n1, n2)

        # det(M) = n2/n1
        @test det(M) ≈ n2/n1
    end

    @testset "Symmetry" begin
        n1, n2 = 1.5, 2.5
        M = PhoXonic.interface_matrix(n1, n2)

        # Matrix is symmetric
        @test M[1, 2] ≈ M[2, 1]
    end
end

@testset "System matrix" begin
    @testset "Single layer" begin
        # Air | Glass (d=100nm) | Air at λ=600nm
        air = Dielectric(1.0)
        glass = Dielectric(2.25)  # n = 1.5
        d = 100.0
        λ = 600.0

        layer = Layer(glass, d)
        ml = Multilayer([layer], air, air)

        M = PhoXonic.system_matrix(ml, λ)

        # Verify structure: M = I_12 * P_2 * I_21
        n_inc = sqrt(1.0)
        n_glass = sqrt(2.25)
        n_sub = sqrt(1.0)

        I12 = PhoXonic.interface_matrix(n_inc, n_glass)
        P2 = PhoXonic.propagation_matrix(n_glass, d, λ)
        I21 = PhoXonic.interface_matrix(n_glass, n_sub)

        M_expected = I12 * P2 * I21
        @test M ≈ M_expected
    end

    @testset "Two layers" begin
        # Air | Material A | Material B | Air
        air = Dielectric(1.0)
        mat_a = Dielectric(4.0)  # n = 2
        mat_b = Dielectric(2.25)  # n = 1.5
        d_a, d_b = 50.0, 75.0
        λ = 500.0

        layers = [Layer(mat_a, d_a), Layer(mat_b, d_b)]
        ml = Multilayer(layers, air, air)

        M = PhoXonic.system_matrix(ml, λ)

        # Check matrix properties
        @test size(M) == (2, 2)
        @test !iszero(det(M))  # Should be non-singular
    end

    @testset "Matrix product order" begin
        # Verify correct ordering: from incident to substrate
        air = Dielectric(1.0)
        glass = Dielectric(2.25)
        d = 100.0
        λ = 600.0

        layer = Layer(glass, d)
        ml = Multilayer([layer], air, air)

        M = PhoXonic.system_matrix(ml, λ)

        # For a single layer, the system matrix relates:
        # [E_inc^+; E_inc^-] = M * [E_sub^+; E_sub^-]
        # where E^+ is forward and E^- is backward traveling
        # M11 should be non-zero (forward transmission exists)
        @test abs(M[1, 1]) > 0
    end
end
