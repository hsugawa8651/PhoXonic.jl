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

    @testset "an unsupported dimension is reported as such" begin
        s3 = Solver(TransverseEM(), geo3, (8, 8, 8), DenseMethod(); cutoff=2)
        ω_values = [0.5, 1.0, 1.5]
        k_points = [[0.1, 0.1, 0.1]]

        # The failure must name the dimension, not blame the algorithm.
        @test_throws ArgumentError compute_dos(s3, ω_values, k_points, DirectGF())
        @test_throws "2D solvers only" compute_dos(s3, ω_values, k_points, DirectGF())
        @test_throws "2D solvers only" compute_dos(s3, ω_values, k_points, RSKGF())
        @test_throws "MatrixFreeGF() only" compute_ldos(
            s3, [0.5, 0.5, 0.5], ω_values, k_points, DirectGF()
        )

        # 1D DOS exists, but only without an explicit GFMethod.
        s1 = Solver(Photonic1D(), geo1, (64,); cutoff=5)
        @test_throws "2D solvers only" compute_dos(s1, ω_values, [0.1, 0.2], DirectGF())
    end

    @testset "a wave vector is the same wave vector at every entry point" begin
        # solve and solve_at_k used to accept different shapes of k, in opposite
        # directions, and a shape that matched neither was reported as an unsupported
        # solver method.  Every shape now reaches every entry point.
        #
        # LOBPCGMethod appears in 2D only.  The 1D fixture here has a matrix dimension
        # of 7, and asking for two bands puts LOBPCG's search block at 3*2 = 6 columns
        # against those 7 rows.  There it breaks down: it throws, or it returns an
        # eigenvalue that is not in the spectrum and reports convergence (#95).

        function failure_message(f)
            try
                f()
            catch e
                return sprint(showerror, e)
            end
            return ""
        end

        @testset "1D: a scalar, a vector and a tuple agree" begin
            for method in (DenseMethod(), KrylovKitMethod())
                s = Solver(Photonic1D(), geo1, (16,), method; cutoff=3)
                scalar, _ = solve(s, 0.3; bands=1:2)

                @test solve(s, [0.3]; bands=1:2)[1] ≈ scalar
                @test solve(s, (0.3,); bands=1:2)[1] ≈ scalar

                # solve_at_k took a vector, but not a scalar with KrylovKitMethod.
                @test solve_at_k(s, 0.3, method; bands=1:2) ≈ scalar
                @test solve_at_k(s, [0.3], method; bands=1:2) ≈ scalar
                @test solve_at_k(s, (0.3,), method; bands=1:2) ≈ scalar
            end
        end

        @testset "2D: an integer Γ, a tuple and a vector agree" begin
            # Normalization is a type conversion, so a deterministic method returns the
            # same numbers whichever shape the caller wrote.  The iterative methods
            # start from random vectors and do not repeat themselves at Γ, where the
            # lowest TM band sits at zero; for those the point is that the call is
            # accepted at all.
            s_dense = Solver(TMWave(), geo2, (16, 16), DenseMethod(); cutoff=3)

            gamma, _ = solve(s_dense, [0.0, 0.0]; bands=1:2)
            @test solve(s_dense, [0, 0]; bands=1:2)[1] == gamma
            @test solve_at_k(s_dense, [0, 0], DenseMethod(); bands=1:2) == gamma

            k, _ = solve(s_dense, [0.1, 0.2]; bands=1:2)
            @test solve(s_dense, (0.1, 0.2); bands=1:2)[1] == k
            @test solve_at_k(s_dense, (0.1, 0.2), DenseMethod(); bands=1:2) == k

            # _solve_krylovkit asks for a Vector{Float64}.  An integer Γ and a tuple
            # used to reach it unconverted, or to reach build_matrices, which has no
            # method for a tuple.
            for method in (KrylovKitMethod(), LOBPCGMethod())
                s = Solver(TMWave(), geo2, (16, 16), method; cutoff=3)

                @test length(solve(s, [0, 0]; bands=1:2)[1]) == 2
                @test length(solve(s, (0.1, 0.2); bands=1:2)[1]) == 2
                @test length(solve_at_k(s, [0, 0], method; bands=1:2)) == 2
                @test length(solve_at_k(s, (0.1, 0.2), method; bands=1:2)) == 2
            end
        end

        @testset "3D: a tuple and a vector agree" begin
            # DenseMethod only.  KrylovKitMethod is the matrix-free path, and it has no
            # operator for TransverseEM, so it throws from inside KrylovKit before the
            # wave vector matters.  LOBPCGMethod starts from random vectors and does not
            # repeat itself, which is what the 2D testset above compares against.
            s = Solver(TransverseEM(), geo3, (8, 8, 8), DenseMethod(); cutoff=2)
            k, _ = solve(s, [0.1, 0.2, 0.3]; bands=1:2)

            @test solve(s, (0.1, 0.2, 0.3); bands=1:2)[1] ≈ k
            @test solve_at_k(s, (0.1, 0.2, 0.3), DenseMethod(); bands=1:2) ≈ k
            @test solve_at_k_with_vectors(s, (0.1, 0.2, 0.3), DenseMethod(); bands=1:2)[1] ≈
                k
        end

        @testset "the wrong number of components names the dimension, not the method" begin
            s1 = Solver(Photonic1D(), geo1, (16,); cutoff=3)
            s2 = Solver(TMWave(), geo2, (16, 16); cutoff=3)
            s3 = Solver(TransverseEM(), geo3, (8, 8, 8), DenseMethod(); cutoff=2)

            @test_throws ArgumentError solve(s2, [0.1]; bands=1:1)
            @test_throws "must have 2 components" solve(s2, [0.1]; bands=1:1)
            @test_throws "must have 2 components" solve(s2, 0.3; bands=1:1)
            @test_throws "must have 3 components" solve_at_k(
                s3, (0.1, 0.2), DenseMethod(); bands=1:1
            )
            @test_throws "must be a real number or a 1-element container" solve(
                s1, [0.1, 0.2]; bands=1:1
            )

            # The message used to blame the solver method, and then list that same
            # method among the supported ones.  It must not do that again.
            for f in (
                () -> solve(s2, [0.1]; bands=1:1),
                () -> solve(s1, [0.1, 0.2]; bands=1:1),
                () -> solve_at_k(s3, (0.1, 0.2), DenseMethod(); bands=1:1),
            )
                @test !occursin("Unsupported solver method", failure_message(f))
            end
        end
    end
end
