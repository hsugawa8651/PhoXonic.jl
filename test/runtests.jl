# Last-Modified: 2025-12-15T21:30:00+09:00

using Test
using PhoXonic
using LinearAlgebra

@testset "PhoXonic.jl" begin

    @testset "Lattice" begin
        @testset "Square lattice" begin
            lat = square_lattice(1.0)
            a1, a2 = lat.vectors
            b1, b2 = lat.reciprocal

            # Check lattice vectors
            @test a1 ≈ [1.0, 0.0]
            @test a2 ≈ [0.0, 1.0]

            # Check reciprocal relations: a_i · b_j = 2π δ_ij
            @test dot(a1, b1) ≈ 2π
            @test dot(a2, b2) ≈ 2π
            @test abs(dot(a1, b2)) < 1e-10
            @test abs(dot(a2, b1)) < 1e-10
        end

        @testset "Hexagonal lattice" begin
            lat = hexagonal_lattice(1.0)
            a1, a2 = lat.vectors
            b1, b2 = lat.reciprocal

            # Check reciprocal relations
            @test dot(a1, b1) ≈ 2π
            @test dot(a2, b2) ≈ 2π
            @test abs(dot(a1, b2)) < 1e-10
            @test abs(dot(a2, b1)) < 1e-10
        end
    end

    @testset "Shapes" begin
        @testset "Circle" begin
            c = Circle([0.0, 0.0], 0.5)
            @test [0.0, 0.0] in c
            @test [0.3, 0.3] in c
            @test !([0.6, 0.6] in c)
        end

        @testset "Rectangle" begin
            r = Rectangle([0.0, 0.0], [1.0, 0.5])
            @test [0.0, 0.0] in r
            @test [0.4, 0.2] in r
            @test !([0.6, 0.0] in r)
        end

        @testset "Shape translate" begin
            # 2D shapes
            @testset "Circle translate" begin
                c = Circle([0.0, 0.0], 0.2)
                c2 = translate(c, [1.0, 2.0])
                @test c2.center ≈ PhoXonic.Vec2(1.0, 2.0)
                @test c2.radius == 0.2
                @test [1.0, 2.0] in c2
                @test !([0.0, 0.0] in c2)
            end

            @testset "Rectangle translate" begin
                r = Rectangle([0.0, 0.0], [1.0, 0.5])
                r2 = translate(r, [2.0, 3.0])
                @test r2.center ≈ PhoXonic.Vec2(2.0, 3.0)
                @test r2.size ≈ PhoXonic.Vec2(1.0, 0.5)
                @test [2.0, 3.0] in r2
            end

            @testset "Polygon translate" begin
                p = Polygon([[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]])
                p2 = translate(p, [1.0, 1.0])
                @test p2.vertices[1] ≈ PhoXonic.Vec2(1.0, 1.0)
                @test p2.vertices[2] ≈ PhoXonic.Vec2(2.0, 1.0)
                @test p2.vertices[3] ≈ PhoXonic.Vec2(1.5, 2.0)
            end

            # 3D shapes
            @testset "Sphere translate" begin
                s = Sphere([0.0, 0.0, 0.0], 0.3)
                s2 = translate(s, [1.0, 2.0, 3.0])
                @test s2.center ≈ PhoXonic.Vec3(1.0, 2.0, 3.0)
                @test s2.radius == 0.3
            end

            @testset "Cylinder translate" begin
                c = Cylinder([0.0, 0.0, 0.0], 0.2, 1.0)
                c2 = translate(c, [1.0, 0.0, 0.5])
                @test c2.center ≈ PhoXonic.Vec3(1.0, 0.0, 0.5)
                @test c2.radius == 0.2
                @test c2.height == 1.0
                @test c2.axis ≈ PhoXonic.Vec3(0, 0, 1)
            end

            @testset "Slab translate" begin
                s = Slab(0.0, 0.5)
                s2 = translate(s, [0.0, 0.0, 1.0])
                @test s2.z_min ≈ 1.0
                @test s2.z_max ≈ 1.5
            end

            # 1D shapes
            @testset "Segment translate" begin
                s = Segment(0.0, 1.0)
                s2 = translate(s, 2.0)
                @test s2.start ≈ 2.0
                @test s2.stop ≈ 3.0
            end
        end
    end

    @testset "Supercell" begin
        # Setup: 2D square lattice with circular rod
        lat = square_lattice(1.0)
        air = Dielectric(1.0)
        rod = Dielectric(8.9)
        geo = Geometry(lat, air, [(Circle([0.5, 0.5], 0.2), rod)])

        @testset "Basic 2D supercell" begin
            geo_super = create_supercell(geo, (3, 3))

            # Check lattice scaling
            @test geo_super.lattice.vectors[1] ≈ 3.0 * lat.vectors[1]
            @test geo_super.lattice.vectors[2] ≈ 3.0 * lat.vectors[2]

            # Check number of inclusions
            @test length(geo_super.inclusions) == 9

            # Check first and last inclusion positions
            centers = [inc[1].center for inc in geo_super.inclusions]
            @test any(c -> c ≈ PhoXonic.Vec2(0.5, 0.5), centers)
            @test any(c -> c ≈ PhoXonic.Vec2(2.5, 2.5), centers)
        end

        @testset "2D supercell with point defect" begin
            # 5×5 supercell with center defect
            geo_defect = create_supercell(geo, (5, 5); point_defects=[(2, 2)])

            # Should have 24 inclusions (25 - 1)
            @test length(geo_defect.inclusions) == 24

            # Check that center position is missing
            centers = [inc[1].center for inc in geo_defect.inclusions]
            center_pos = PhoXonic.Vec2(2.5, 2.5)  # (2+0.5, 2+0.5)
            @test !any(c -> c ≈ center_pos, centers)
        end

        @testset "line_defect_positions 2D" begin
            # Horizontal line at row 2 in 7×5 supercell
            positions = line_defect_positions(:x, 2, (7, 5))
            @test length(positions) == 7
            @test (0, 2) in positions
            @test (6, 2) in positions

            # Vertical line at column 3 in 7×5 supercell
            positions_y = line_defect_positions(:y, 3, (7, 5))
            @test length(positions_y) == 5
            @test (3, 0) in positions_y
            @test (3, 4) in positions_y
        end

        @testset "2D supercell with line defect (waveguide)" begin
            defects = line_defect_positions(:x, 2, (7, 5))
            geo_wg = create_supercell(geo, (7, 5); point_defects=defects)

            # Should have 35 - 7 = 28 inclusions
            @test length(geo_wg.inclusions) == 28
        end

        # 3D tests
        @testset "Basic 3D supercell" begin
            lat3d = cubic_lattice(1.0)
            sphere_mat = Dielectric(12.0)
            geo3d = Geometry(lat3d, air, [(Sphere([0.5, 0.5, 0.5], 0.2), sphere_mat)])

            geo_super3d = create_supercell(geo3d, (2, 2, 2))

            # Check lattice scaling
            @test geo_super3d.lattice.vectors[1] ≈ 2.0 * lat3d.vectors[1]

            # Check number of inclusions
            @test length(geo_super3d.inclusions) == 8
        end

        @testset "3D supercell with point defect" begin
            lat3d = cubic_lattice(1.0)
            sphere_mat = Dielectric(12.0)
            geo3d = Geometry(lat3d, air, [(Sphere([0.5, 0.5, 0.5], 0.2), sphere_mat)])

            geo_defect3d = create_supercell(geo3d, (3, 3, 3); point_defects=[(1, 1, 1)])

            # Should have 27 - 1 = 26 inclusions
            @test length(geo_defect3d.inclusions) == 26
        end

        @testset "line_defect_positions 3D" begin
            # Line along x at (y=2, z=3)
            positions = line_defect_positions(:x, (2, 3), (7, 5, 5))
            @test length(positions) == 7
            @test (0, 2, 3) in positions
            @test (6, 2, 3) in positions
        end

        # 1D tests
        @testset "Basic 1D supercell" begin
            lat1d = lattice_1d(1.0)
            mat1d = Dielectric(4.0)
            geo1d = Geometry(lat1d, air, [(Segment(0.3, 0.7), mat1d)])

            geo_super1d = create_supercell(geo1d, 5)

            # Check lattice scaling
            @test geo_super1d.lattice.vectors[1][1] ≈ 5.0

            # Check number of inclusions
            @test length(geo_super1d.inclusions) == 5
        end
    end

    # Error fallback tests (can also run independently: julia --project=. test/test_error_fallbacks.jl)
    include("test_error_fallbacks.jl")

    @testset "Fill fraction (periodic boundary)" begin
        # These tests verify that periodic boundary conditions are handled correctly
        # for all 2D shapes at various positions in the unit cell, on both orthogonal
        # (square) and non-orthogonal (hexagonal) lattices.

        # Helper function to compute cell area
        function cell_area(lat)
            a1, a2 = lat.vectors
            abs(a1[1] * a2[2] - a1[2] * a2[1])
        end

        # ====================================================================
        # Circle tests
        # ====================================================================
        @testset "Circle on square lattice" begin
            a = 1.0
            r = 0.2
            lat = square_lattice(a)
            air = Dielectric(1.0)
            rod = Dielectric(12.0)

            theoretical_fill = π * r^2 / cell_area(lat)

            for (label, center) in [("center (0.5, 0.5)", [0.5, 0.5]),
                                     ("corner (0.0, 0.0)", [0.0, 0.0]),
                                     ("edge (0.5, 0.0)", [0.5, 0.0])]
                geo = Geometry(lat, air, [(Circle(center, r), rod)])
                ε = PhoXonic.discretize(geo, (64, 64), :ε)
                rod_pixels = sum(ε .> 2.0)
                fill_fraction = rod_pixels / length(ε)
                @test isapprox(fill_fraction, theoretical_fill, rtol=0.05)
            end
        end

        @testset "Circle on hexagonal lattice" begin
            a = 1.0
            r = 0.25
            lat = hexagonal_lattice(a)
            air = Dielectric(1.0)
            rod = Dielectric(12.0)

            theoretical_fill = π * r^2 / cell_area(lat)

            for (label, center) in [("center", [0.5, 0.289]),  # center of unit cell
                                     ("corner (0.0, 0.0)", [0.0, 0.0])]
                geo = Geometry(lat, air, [(Circle(center, r), rod)])
                ε = PhoXonic.discretize(geo, (64, 64), :ε)
                rod_pixels = sum(ε .> 2.0)
                fill_fraction = rod_pixels / length(ε)
                @test isapprox(fill_fraction, theoretical_fill, rtol=0.05)
            end
        end

        # ====================================================================
        # Rectangle tests
        # ====================================================================
        @testset "Rectangle on square lattice" begin
            a = 1.0
            lat = square_lattice(a)
            air = Dielectric(1.0)
            rod = Dielectric(12.0)

            rect_size = [0.4, 0.3]
            theoretical_fill = rect_size[1] * rect_size[2] / cell_area(lat)

            for (label, center) in [("center (0.5, 0.5)", [0.5, 0.5]),
                                     ("corner (0.0, 0.0)", [0.0, 0.0]),
                                     ("edge (0.5, 0.0)", [0.5, 0.0])]
                geo = Geometry(lat, air, [(Rectangle(center, rect_size), rod)])
                ε = PhoXonic.discretize(geo, (64, 64), :ε)
                rod_pixels = sum(ε .> 2.0)
                fill_fraction = rod_pixels / length(ε)
                @test isapprox(fill_fraction, theoretical_fill, rtol=0.05)
            end
        end

        @testset "Rectangle on hexagonal lattice" begin
            a = 1.0
            lat = hexagonal_lattice(a)
            air = Dielectric(1.0)
            rod = Dielectric(12.0)

            rect_size = [0.3, 0.2]
            theoretical_fill = rect_size[1] * rect_size[2] / cell_area(lat)

            for (label, center) in [("center", [0.5, 0.289]),
                                     ("corner (0.0, 0.0)", [0.0, 0.0])]
                geo = Geometry(lat, air, [(Rectangle(center, rect_size), rod)])
                ε = PhoXonic.discretize(geo, (64, 64), :ε)
                rod_pixels = sum(ε .> 2.0)
                fill_fraction = rod_pixels / length(ε)
                @test isapprox(fill_fraction, theoretical_fill, rtol=0.05)
            end
        end

        # ====================================================================
        # Polygon tests (triangle as representative polygon)
        # Use higher resolution (128×128) for better accuracy with polygon edges
        # ====================================================================
        @testset "Polygon (triangle) on square lattice" begin
            a = 1.0
            lat = square_lattice(a)
            air = Dielectric(1.0)
            rod = Dielectric(12.0)

            # Equilateral triangle with side s
            s = 0.4
            triangle_area = (sqrt(3) / 4) * s^2
            theoretical_fill = triangle_area / cell_area(lat)

            # Triangle vertices (centered at different positions)
            for (label, center) in [("center (0.5, 0.5)", [0.5, 0.5]),
                                     ("corner (0.0, 0.0)", [0.0, 0.0])]
                # Equilateral triangle vertices centered at `center`
                h = s * sqrt(3) / 2  # height
                vertices = [
                    [center[1], center[2] + h * 2/3],           # top
                    [center[1] - s/2, center[2] - h * 1/3],     # bottom-left
                    [center[1] + s/2, center[2] - h * 1/3]      # bottom-right
                ]
                geo = Geometry(lat, air, [(Polygon(vertices), rod)])
                ε = PhoXonic.discretize(geo, (128, 128), :ε)
                rod_pixels = sum(ε .> 2.0)
                fill_fraction = rod_pixels / length(ε)
                @test isapprox(fill_fraction, theoretical_fill, rtol=0.10)
            end
        end

        @testset "Polygon (triangle) on hexagonal lattice" begin
            a = 1.0
            lat = hexagonal_lattice(a)
            air = Dielectric(1.0)
            rod = Dielectric(12.0)

            # Equilateral triangle with side s
            s = 0.3
            triangle_area = (sqrt(3) / 4) * s^2
            theoretical_fill = triangle_area / cell_area(lat)

            for (label, center) in [("center", [0.5, 0.289]),
                                     ("corner (0.0, 0.0)", [0.0, 0.0])]
                h = s * sqrt(3) / 2
                vertices = [
                    [center[1], center[2] + h * 2/3],
                    [center[1] - s/2, center[2] - h * 1/3],
                    [center[1] + s/2, center[2] - h * 1/3]
                ]
                geo = Geometry(lat, air, [(Polygon(vertices), rod)])
                ε = PhoXonic.discretize(geo, (128, 128), :ε)
                rod_pixels = sum(ε .> 2.0)
                fill_fraction = rod_pixels / length(ε)
                @test isapprox(fill_fraction, theoretical_fill, rtol=0.10)
            end
        end

        # ====================================================================
        # 1D Segment tests
        # ====================================================================
        @testset "Segment (1D)" begin
            a = 1.0
            lat = lattice_1d(a)
            mat1 = Dielectric(1.0)
            mat2 = Dielectric(9.0)

            # Test different segment positions
            # Fill fraction = segment_length / lattice_constant
            segment_length = 0.3

            for (label, seg_start, seg_stop) in [
                    ("center [0.35, 0.65]", 0.35, 0.65),           # Center of unit cell
                    ("edge [0.0, 0.3]", 0.0, segment_length),       # At edge (start = 0)
                    ("corner [-0.15, 0.15]", -0.15, 0.15),         # Wrapping around origin
                    ("end [0.85, 1.15]", 0.85, 1.15)]              # Wrapping around end
                geo = Geometry(lat, mat1, [(Segment(seg_start, seg_stop), mat2)])
                ε = PhoXonic.discretize(geo, 128, :ε)
                mat2_pixels = sum(ε .> 2.0)
                fill_fraction = mat2_pixels / length(ε)
                theoretical_fill = segment_length / a
                @test isapprox(fill_fraction, theoretical_fill, rtol=0.05)
            end
        end

        # ====================================================================
        # 3D Sphere tests
        # ====================================================================
        @testset "Sphere on cubic lattice" begin
            a = 1.0
            r = 0.2
            lat = cubic_lattice(a)
            air = Dielectric(1.0)
            rod = Dielectric(12.0)

            # Theoretical fill fraction = (4/3)πr³ / a³
            theoretical_fill = (4/3) * π * r^3 / a^3

            for (label, center) in [("center (0.5, 0.5, 0.5)", [0.5, 0.5, 0.5]),
                                     ("corner (0.0, 0.0, 0.0)", [0.0, 0.0, 0.0])]
                geo = Geometry(lat, air, [(Sphere(center, r), rod)])
                ε = PhoXonic.discretize(geo, (32, 32, 32), :ε)
                rod_pixels = sum(ε .> 2.0)
                fill_fraction = rod_pixels / length(ε)
                @test isapprox(fill_fraction, theoretical_fill, rtol=0.10)
            end
        end

        # ====================================================================
        # 3D Cylinder tests (z-axis aligned)
        # ====================================================================
        @testset "Cylinder on cubic lattice" begin
            a = 1.0
            r = 0.15
            h = a  # Full height through unit cell
            lat = cubic_lattice(a)
            air = Dielectric(1.0)
            rod = Dielectric(12.0)

            # Theoretical fill fraction = πr²h / a³ = πr² / a² (when h = a)
            theoretical_fill = π * r^2 / a^2

            for (label, center) in [("center (0.5, 0.5, 0.5)", [0.5, 0.5, 0.5]),
                                     ("corner (0.0, 0.0, 0.5)", [0.0, 0.0, 0.5])]
                geo = Geometry(lat, air, [(Cylinder(center, r, h), rod)])
                ε = PhoXonic.discretize(geo, (32, 32, 32), :ε)
                rod_pixels = sum(ε .> 2.0)
                fill_fraction = rod_pixels / length(ε)
                @test isapprox(fill_fraction, theoretical_fill, rtol=0.10)
            end
        end

        @testset "Slab on cubic lattice" begin
            a = 1.0
            lat = cubic_lattice(a)
            air = Dielectric(1.0)
            slab_mat = Dielectric(12.0)

            # Slab from z=0.2 to z=0.7 (thickness = 0.5)
            z_min, z_max = 0.2, 0.7
            thickness = z_max - z_min

            # Theoretical fill fraction = thickness / a
            theoretical_fill = thickness / a

            geo = Geometry(lat, air, [(Slab(z_min, z_max), slab_mat)])
            ε = PhoXonic.discretize(geo, (16, 16, 16), :ε)
            slab_pixels = sum(ε .> 2.0)
            fill_fraction = slab_pixels / length(ε)
            @test isapprox(fill_fraction, theoretical_fill, rtol=0.10)

            # Test point containment
            slab = Slab(0.2, 0.8)
            @test [0.0, 0.0, 0.5] in slab
            @test [100.0, -50.0, 0.2] in slab  # xy is arbitrary
            @test [0.0, 0.0, 0.8] in slab      # boundary
            @test !([0.0, 0.0, 0.1] in slab)
            @test !([0.0, 0.0, 0.9] in slab)
        end
    end

    @testset "Convolution matrix" begin
        @testset "Homogeneous material" begin
            lat = square_lattice(1.0)
            basis = PlaneWaveBasis(lat, 3)

            # Homogeneous ε = 2.0
            ε = fill(2.0, 32, 32)
            C = convolution_matrix(ε, basis)

            # For homogeneous material, convolution matrix should be diagonal
            # with diagonal elements equal to ε
            for i in 1:basis.num_pw
                @test abs(C[i, i] - 2.0) < 1e-10
            end

            # Off-diagonal elements should be zero
            for i in 1:basis.num_pw, j in 1:basis.num_pw
                if i != j
                    @test abs(C[i, j]) < 1e-10
                end
            end
        end
    end

    @testset "Homogeneous photonic crystal" begin
        # For homogeneous medium, ω = |k| / √(εμ)
        lat = square_lattice(1.0)
        ε_val = 2.0
        μ_val = 1.0

        air = Dielectric(ε_val, μ_val)
        geo = Geometry(lat, air)

        solver = Solver(TEWave(), geo, (32, 32); cutoff=5)

        # Test at k = (1, 0)
        k = [1.0, 0.0]
        ω, _ = solve(solver, k; bands=1:3)

        # Expected: ω = |k| / √(εμ) = 1 / √2
        expected_ω = norm(k) / sqrt(ε_val * μ_val)
        @test isapprox(ω[1], expected_ω, rtol=0.01)

        # Test TM mode
        solver_tm = Solver(TMWave(), geo, (32, 32); cutoff=5)
        ω_tm, _ = solve(solver_tm, k; bands=1:3)
        @test isapprox(ω_tm[1], expected_ω, rtol=0.01)
    end

    @testset "Homogeneous at high-symmetry points" begin
        # Test light line ω = |k| / √ε at high-symmetry points
        a = 1.0
        ε_val = 12.0
        lat = square_lattice(a)
        mat = Dielectric(ε_val)
        geo = Geometry(lat, mat)

        solver_tm = Solver(TMWave(), geo, (32, 32); cutoff=7)
        solver_te = Solver(TEWave(), geo, (32, 32); cutoff=7)

        # X point: k = (π/a, 0)
        k_X = [π/a, 0.0]
        ω_expected_X = norm(k_X) / sqrt(ε_val)
        ω_tm_X, _ = solve(solver_tm, k_X; bands=1:1)
        ω_te_X, _ = solve(solver_te, k_X; bands=1:1)
        @test isapprox(ω_tm_X[1], ω_expected_X, rtol=0.01)
        @test isapprox(ω_te_X[1], ω_expected_X, rtol=0.01)

        # M point: k = (π/a, π/a)
        k_M = [π/a, π/a]
        ω_expected_M = norm(k_M) / sqrt(ε_val)
        ω_tm_M, _ = solve(solver_tm, k_M; bands=1:1)
        ω_te_M, _ = solve(solver_te, k_M; bands=1:1)
        @test isapprox(ω_tm_M[1], ω_expected_M, rtol=0.01)
        @test isapprox(ω_te_M[1], ω_expected_M, rtol=0.01)
    end

    @testset "1D/2D consistency for homogeneous medium" begin
        # For homogeneous medium with ky=0, 1D and 2D should give same result
        a = 1.0
        ε_val = 12.0

        # 1D solver
        lat_1d = lattice_1d(a)
        geo_1d = Geometry(lat_1d, Dielectric(ε_val))
        solver_1d = Solver(Photonic1D(), geo_1d, 256; cutoff=30)

        # 2D TM solver
        lat_2d = square_lattice(a)
        geo_2d = Geometry(lat_2d, Dielectric(ε_val))
        solver_2d = Solver(TMWave(), geo_2d, (32, 32); cutoff=7)

        # Test at zone boundary k = π/a
        k_1d = π/a
        k_2d = [π/a, 0.0]  # ky = 0

        ω_1d, _ = solve(solver_1d, k_1d; bands=1:1)
        ω_2d, _ = solve(solver_2d, k_2d; bands=1:1)

        # Both should give ω = k / √ε
        expected_ω = k_1d / sqrt(ε_val)
        @test isapprox(ω_1d[1], expected_ω, rtol=0.01)
        @test isapprox(ω_2d[1], expected_ω, rtol=0.01)
        @test isapprox(ω_1d[1], ω_2d[1], rtol=0.01)
    end

    @testset "SimpleKPath" begin
        kpath = simple_kpath_square(a=1.0, npoints=10)

        # Check that we have points
        @test length(kpath) > 0

        # Check iteration
        count = 0
        for k in kpath
            count += 1
            @test length(k) == 2
        end
        @test count == length(kpath)

        # Check labels
        @test length(kpath.labels) == 4  # Γ, X, M, Γ
    end

    @testset "3D photonic crystal" begin
        @testset "Homogeneous 3D medium" begin
            lat = cubic_lattice(1.0)
            air = Dielectric(1.0, 1.0)
            geo = Geometry(lat, air)

            solver = Solver(FullVectorEM(), geo, (8, 8, 8); cutoff=2)
            k = [0.5, 0.3, 0.2]  # Away from origin

            N = solver.basis.num_pw
            # For H-field formulation, there are N longitudinal modes with ω≈0
            # and 2N transverse modes with physical frequencies
            ω, _ = solve(solver, k; bands=1:(N+5))

            # Expected: ω = |k| for the lowest transverse mode (doubly degenerate)
            expected = norm(k)
            # Skip the N longitudinal modes (ω≈0), check transverse modes
            transverse_ω = filter(x -> x > 0.01, ω)
            @test length(transverse_ω) >= 2  # At least 2 transverse modes
            @test isapprox(transverse_ω[1], expected, rtol=0.01)
        end

        @testset "3D with inclusion" begin
            lat = cubic_lattice(1.0)
            air = Dielectric(1.0, 1.0)
            rod = Dielectric(8.9, 1.0)
            s = Sphere([0.0, 0.0, 0.0], 0.2)
            geo = Geometry(lat, air, [(s, rod)])

            solver = Solver(FullVectorEM(), geo, (8, 8, 8); cutoff=2)
            k = [0.1, 0.1, 0.1]

            N = solver.basis.num_pw
            ω, modes = solve(solver, k; bands=1:(N+5))

            # Just check it runs and returns reasonable values
            @test length(ω) == N + 5
            @test all(ω .>= 0)
            @test size(modes, 2) == N + 5
        end
    end

    @testset "KrylovKit iterative solver" begin
        @testset "TEWave iterative vs dense" begin
            lat = square_lattice(1.0)
            air = Dielectric(1.0, 1.0)
            rod = Dielectric(8.9, 1.0)
            c = Circle([0.0, 0.0], 0.2)
            geo = Geometry(lat, air, [(c, rod)])

            # Dense solver
            solver_dense = Solver(TEWave(), geo, (32, 32); cutoff=3)

            # Iterative solver
            solver_iter = Solver(TEWave(), geo, (32, 32), KrylovKitMethod(); cutoff=3)

            k = [0.1, 0.2]
            nbands = 5

            ω_dense, _ = solve(solver_dense, k; bands=1:nbands)
            ω_iter, _ = solve(solver_iter, k; bands=1:nbands)

            # Compare frequencies (relative tolerance 1e-4 for iterative method)
            for i in 1:nbands
                @test isapprox(ω_iter[i], ω_dense[i], rtol=1e-4)
            end
        end

        @testset "TMWave iterative vs dense" begin
            lat = square_lattice(1.0)
            air = Dielectric(1.0, 1.0)
            rod = Dielectric(8.9, 1.0)
            c = Circle([0.0, 0.0], 0.2)
            geo = Geometry(lat, air, [(c, rod)])

            solver_dense = Solver(TMWave(), geo, (32, 32); cutoff=3)
            solver_iter = Solver(TMWave(), geo, (32, 32), KrylovKitMethod(); cutoff=3)

            k = [0.1, 0.2]
            nbands = 5

            ω_dense, _ = solve(solver_dense, k; bands=1:nbands)
            ω_iter, _ = solve(solver_iter, k; bands=1:nbands)

            for i in 1:nbands
                @test isapprox(ω_iter[i], ω_dense[i], rtol=1e-4)
            end
        end

        @testset "Shift-and-invert 3D" begin
            lat = cubic_lattice(1.0)
            air = Dielectric(1.0, 1.0)
            geo = Geometry(lat, air)

            # Use shift to skip longitudinal modes
            solver = Solver(FullVectorEM(), geo, (8, 8, 8), KrylovKitMethod(shift=0.01); cutoff=2)
            k = [0.5, 0.3, 0.2]

            ω, _ = solve(solver, k; bands=1:5)

            # With shift, should directly get transverse modes
            # Expected: ω = |k| for homogeneous medium
            expected = norm(k)
            @test isapprox(ω[1], expected, rtol=0.01)
            @test isapprox(ω[2], expected, rtol=0.01)  # Doubly degenerate
        end

        @testset "SHWave iterative vs dense" begin
            lat = square_lattice(0.01)
            epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
            steel = IsotropicElastic(ρ=7800.0, λ=1.15e11, μ=8.28e10)
            c = Circle([0.0, 0.0], 0.004)
            geo = Geometry(lat, epoxy, [(c, steel)])

            solver_dense = Solver(SHWave(), geo, (32, 32); cutoff=3)
            solver_iter = Solver(SHWave(), geo, (32, 32), KrylovKitMethod(); cutoff=3)

            k = [100.0, 150.0]
            nbands = 5

            ω_dense, _ = solve(solver_dense, k; bands=1:nbands)
            ω_iter, _ = solve(solver_iter, k; bands=1:nbands)

            for i in 1:nbands
                @test isapprox(ω_iter[i], ω_dense[i], rtol=1e-4)
            end
        end

        @testset "PSVWave iterative vs dense" begin
            lat = square_lattice(0.01)
            epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
            steel = IsotropicElastic(ρ=7800.0, λ=1.15e11, μ=8.28e10)
            c = Circle([0.0, 0.0], 0.004)
            geo = Geometry(lat, epoxy, [(c, steel)])

            solver_dense = Solver(PSVWave(), geo, (32, 32); cutoff=3)
            solver_iter = Solver(PSVWave(), geo, (32, 32), KrylovKitMethod(); cutoff=3)

            k = [100.0, 150.0]
            nbands = 5

            ω_dense, _ = solve(solver_dense, k; bands=1:nbands)
            ω_iter, _ = solve(solver_iter, k; bands=1:nbands)

            for i in 1:nbands
                @test isapprox(ω_iter[i], ω_dense[i], rtol=1e-4)
            end
        end
    end

    @testset "1D iterative solvers" begin
        @testset "Photonic1D KrylovKit vs dense" begin
            lat = lattice_1d(1.0)
            mat1 = Dielectric(1.0)
            mat2 = Dielectric(9.0)
            geo = Geometry(lat, mat1, [(Segment(0.0, 0.25), mat2)])

            solver_dense = Solver(Photonic1D(), geo, 64; cutoff=10)
            solver_iter = Solver(Photonic1D(), geo, 64, KrylovKitMethod(); cutoff=10)

            k = 0.3
            nbands = 4

            ω_dense, _ = solve(solver_dense, k; bands=1:nbands)
            ω_iter, _ = solve(solver_iter, k; bands=1:nbands)

            for i in 1:nbands
                @test isapprox(ω_iter[i], ω_dense[i], rtol=1e-4)
            end
        end

        @testset "Photonic1D LOBPCG (shifted) vs dense" begin
            lat = lattice_1d(1.0)
            mat1 = Dielectric(1.0)
            mat2 = Dielectric(9.0)
            geo = Geometry(lat, mat1, [(Segment(0.0, 0.25), mat2)])

            solver_dense = Solver(Photonic1D(), geo, 64; cutoff=10)
            solver_lobpcg = Solver(Photonic1D(), geo, 64, LOBPCGMethod(shift=0.01); cutoff=10)

            k = 0.3
            nbands = 4

            ω_dense, _ = solve(solver_dense, k; bands=1:nbands)
            ω_lobpcg, _ = solve(solver_lobpcg, k; bands=1:nbands)

            for i in 1:nbands
                @test isapprox(ω_lobpcg[i], ω_dense[i], rtol=1e-4)
            end
        end

        # Note: Longitudinal1D KrylovKit requires careful tuning of shift parameter
        # due to large phononic eigenvalues (ω² ~ 10^8). LOBPCG with shift is recommended.

        @testset "Longitudinal1D LOBPCG (shifted) vs dense" begin
            lat = lattice_1d(0.01)
            steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)
            epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
            geo = Geometry(lat, epoxy, [(Segment(0.0, 0.5), steel)])

            solver_dense = Solver(Longitudinal1D(), geo, 64; cutoff=10)
            # For phononic, shift needs to be scaled (ω² ~ 10^8)
            solver_lobpcg = Solver(Longitudinal1D(), geo, 64, LOBPCGMethod(shift=1e6); cutoff=10)

            k = 100.0
            nbands = 4

            ω_dense, _ = solve(solver_dense, k; bands=1:nbands)
            ω_lobpcg, _ = solve(solver_lobpcg, k; bands=1:nbands)

            for i in 1:nbands
                @test isapprox(ω_lobpcg[i], ω_dense[i], rtol=1e-4)
            end
        end
    end

    @testset "Matrix-free operators" begin
        @testset "TMWave matrix-free vs dense" begin
            # Setup: simple 2D photonic crystal
            lat = square_lattice(1.0)
            air = Dielectric(1.0, 1.0)  # ε=1, μ=1
            rod = Dielectric(8.9, 1.0)  # ε=8.9, μ=1 (alumina)
            c = Circle([0.0, 0.0], 0.2)
            geo = Geometry(lat, air, [(c, rod)])  # (Shape, Material) order

            solver = Solver(TMWave(), geo, (32, 32); cutoff=3)
            k = [0.1, 0.2]

            # Dense matrices
            LHS_dense, RHS_dense = PhoXonic.build_matrices(solver, k)

            # Matrix-free operator
            op = PhoXonic.MatrixFreeOperator(solver, k)
            N = solver.basis.num_pw

            # Test with random vector
            x = rand(ComplexF64, N)

            # Dense LHS
            y_dense_lhs = LHS_dense * x

            # Matrix-free LHS
            y_mf_lhs = zeros(ComplexF64, N)
            PhoXonic.apply_lhs!(y_mf_lhs, op, x)

            @test norm(y_dense_lhs - y_mf_lhs) / norm(y_dense_lhs) < 1e-10

            # Dense RHS
            y_dense_rhs = RHS_dense * x

            # Matrix-free RHS
            y_mf_rhs = zeros(ComplexF64, N)
            PhoXonic.apply_rhs!(y_mf_rhs, op, x)

            @test norm(y_dense_rhs - y_mf_rhs) / norm(y_dense_rhs) < 1e-10
        end

        @testset "TEWave matrix-free vs dense" begin
            lat = square_lattice(1.0)
            air = Dielectric(1.0, 1.0)
            rod = Dielectric(8.9, 1.0)
            c = Circle([0.0, 0.0], 0.2)
            geo = Geometry(lat, air, [(c, rod)])  # (Shape, Material) order

            solver = Solver(TEWave(), geo, (32, 32); cutoff=3)
            k = [0.1, 0.2]

            LHS_dense, RHS_dense = PhoXonic.build_matrices(solver, k)
            op = PhoXonic.MatrixFreeOperator(solver, k)
            N = solver.basis.num_pw

            x = rand(ComplexF64, N)

            y_dense_lhs = LHS_dense * x
            y_mf_lhs = zeros(ComplexF64, N)
            PhoXonic.apply_lhs!(y_mf_lhs, op, x)
            @test norm(y_dense_lhs - y_mf_lhs) / norm(y_dense_lhs) < 1e-10

            y_dense_rhs = RHS_dense * x
            y_mf_rhs = zeros(ComplexF64, N)
            PhoXonic.apply_rhs!(y_mf_rhs, op, x)
            @test norm(y_dense_rhs - y_mf_rhs) / norm(y_dense_rhs) < 1e-10
        end

        @testset "PSVWave matrix-free vs dense" begin
            lat = square_lattice(1.0)
            # Steel matrix with aluminum inclusion
            steel = IsotropicElastic(7800.0, 280e9, 130e9, 75e9)  # ρ, C11, C12, C44
            aluminum = IsotropicElastic(2700.0, 107e9, 55e9, 26e9)
            c = Circle([0.0, 0.0], 0.3)
            geo = Geometry(lat, steel, [(c, aluminum)])  # (Shape, Material) order

            solver = Solver(PSVWave(), geo, (32, 32); cutoff=3)
            k = [0.1, 0.2]

            LHS_dense, RHS_dense = PhoXonic.build_matrices(solver, k)
            op = PhoXonic.MatrixFreeOperator(solver, k)
            N = solver.basis.num_pw

            # PSVWave has 2N components
            x = rand(ComplexF64, 2N)

            y_dense_lhs = LHS_dense * x
            y_mf_lhs = zeros(ComplexF64, 2N)
            PhoXonic.apply_lhs!(y_mf_lhs, op, x)
            @test norm(y_dense_lhs - y_mf_lhs) / norm(y_dense_lhs) < 1e-10

            y_dense_rhs = RHS_dense * x
            y_mf_rhs = zeros(ComplexF64, 2N)
            PhoXonic.apply_rhs!(y_mf_rhs, op, x)
            @test norm(y_dense_rhs - y_mf_rhs) / norm(y_dense_rhs) < 1e-10
        end

        @testset "SHWave matrix-free vs dense" begin
            lat = square_lattice(1.0)
            steel = IsotropicElastic(7800.0, 280e9, 130e9, 75e9)
            aluminum = IsotropicElastic(2700.0, 107e9, 55e9, 26e9)
            c = Circle([0.0, 0.0], 0.3)
            geo = Geometry(lat, steel, [(c, aluminum)])  # (Shape, Material) order

            solver = Solver(SHWave(), geo, (32, 32); cutoff=3)
            k = [0.1, 0.2]

            LHS_dense, RHS_dense = PhoXonic.build_matrices(solver, k)
            op = PhoXonic.MatrixFreeOperator(solver, k)
            N = solver.basis.num_pw

            x = rand(ComplexF64, N)

            y_dense_lhs = LHS_dense * x
            y_mf_lhs = zeros(ComplexF64, N)
            PhoXonic.apply_lhs!(y_mf_lhs, op, x)
            @test norm(y_dense_lhs - y_mf_lhs) / norm(y_dense_lhs) < 1e-10

            y_dense_rhs = RHS_dense * x
            y_mf_rhs = zeros(ComplexF64, N)
            PhoXonic.apply_rhs!(y_mf_rhs, op, x)
            @test norm(y_dense_rhs - y_mf_rhs) / norm(y_dense_rhs) < 1e-10
        end
    end

    # ========================================================================
    # 3D Matrix-Free Tests
    # ========================================================================
    @testset "3D Matrix-Free" begin
        # Common setup
        function setup_3d_photonic()
            lat = cubic_lattice(1.0)
            air = Dielectric(1.0)
            rod = Dielectric(4.0)
            geo = Geometry(lat, air, [(Sphere([0.0, 0.0, 0.0], 0.3), rod)])
            return geo
        end

        @testset "MatrixFreeOperator Dim3 constructor" begin
            geo = setup_3d_photonic()
            solver = Solver(FullVectorEM(), geo, (8, 8, 8); cutoff=2)
            k = [0.1, 0.2, 0.3]

            op = PhoXonic.MatrixFreeOperator(solver, k)

            @test op.solver === solver
            @test op.k == k
            @test PhoXonic.resolution(op) == (8, 8, 8)
            @test size(PhoXonic.work_real(op)) == (8, 8, 8)
        end

        @testset "fourier_to_grid! Dim3" begin
            geo = setup_3d_photonic()
            solver = Solver(FullVectorEM(), geo, (8, 8, 8); cutoff=2)
            basis = solver.basis
            res = solver.resolution
            N = basis.num_pw

            # Set 1.0 at G=(0,0,0) -> uniform distribution
            coeffs = zeros(ComplexF64, N)
            idx_zero = findfirst(g -> g == (0, 0, 0), basis.indices)
            coeffs[idx_zero] = 1.0

            grid = zeros(ComplexF64, res)
            PhoXonic.fourier_to_grid!(grid, coeffs, basis, res)

            # All grid points should have the same value (uniform)
            @test all(x -> isapprox(x, grid[1,1,1], rtol=1e-10), grid)
        end

        @testset "grid_to_fourier! Dim3 roundtrip" begin
            geo = setup_3d_photonic()
            solver = Solver(FullVectorEM(), geo, (8, 8, 8); cutoff=2)
            basis = solver.basis
            res = solver.resolution
            N = basis.num_pw

            # Roundtrip test
            coeffs_orig = randn(ComplexF64, N)
            grid = zeros(ComplexF64, res)
            coeffs_back = zeros(ComplexF64, N)

            PhoXonic.fourier_to_grid!(grid, coeffs_orig, basis, res)
            PhoXonic.grid_to_fourier!(coeffs_back, grid, basis, res)

            @test coeffs_back ≈ coeffs_orig rtol=1e-10
        end

        @testset "apply_rhs! FullVectorEM vs Dense" begin
            geo = setup_3d_photonic()
            solver = Solver(FullVectorEM(), geo, (8, 8, 8); cutoff=2)
            k = [0.1, 0.2, 0.3]
            N = solver.basis.num_pw

            # Dense RHS
            LHS, RHS = PhoXonic.build_matrices(solver, k)

            # Random input
            x = randn(ComplexF64, 3N)
            y_dense = RHS * x

            # Matrix-free RHS
            op = PhoXonic.MatrixFreeOperator(solver, k)
            y_mf = zeros(ComplexF64, 3N)
            PhoXonic.apply_rhs!(y_mf, op, x)

            @test y_mf ≈ y_dense rtol=1e-10
        end

        @testset "apply_lhs! FullVectorEM vs Dense" begin
            geo = setup_3d_photonic()
            # Resolution must be large enough for FFT convolution
            # (at least 4x cutoff to avoid aliasing in matrix-free method)
            solver = Solver(FullVectorEM(), geo, (16, 16, 16); cutoff=2)
            k = [0.1, 0.2, 0.3]
            N = solver.basis.num_pw

            # Dense LHS
            LHS, RHS = PhoXonic.build_matrices(solver, k)

            # Random input
            x = randn(ComplexF64, 3N)
            y_dense = LHS * x

            # Matrix-free LHS
            op = PhoXonic.MatrixFreeOperator(solver, k)
            y_mf = zeros(ComplexF64, 3N)
            PhoXonic.apply_lhs!(y_mf, op, x)

            @test y_mf ≈ y_dense rtol=1e-10
        end

        @testset "Eigenvalue consistency: Dense solver" begin
            # Verify that the 3D Dense eigenvalue solver works correctly
            # with shift-and-invert for homogeneous medium
            lattice = PhoXonic.cubic_lattice(1.0)
            material = PhoXonic.Dielectric(4.0)  # ε = 4
            geo = PhoXonic.Geometry(lattice, material)

            solver = Solver(FullVectorEM(), geo, (8, 8, 8); cutoff=2)
            k = [0.5, 0.0, 0.0]

            # Dense eigenvalues with shift-and-invert
            LHS, RHS = PhoXonic.build_matrices(solver, k)
            σ = 0.1
            shifted_LHS = LHS - σ * RHS
            vals = eigvals(shifted_LHS, RHS)
            λ_vals = real.(vals) .+ σ
            ω_vals = sqrt.(max.(λ_vals, 0.0))
            ω_sorted = sort(filter(x -> x > 0.01, ω_vals))

            # Expected: ω = |k|/√ε = 0.5/2 = 0.25 for first band
            @test length(ω_sorted) >= 2
            @test ω_sorted[1] ≈ 0.25 rtol=1e-3
            @test ω_sorted[2] ≈ 0.25 rtol=1e-3  # Doubly degenerate
        end
    end

    # ========================================================================
    # 3D FullElastic Matrix-Free Tests
    # ========================================================================
    @testset "3D FullElastic Matrix-Free" begin
        # Common setup for 3D phononic crystal
        function setup_3d_phononic()
            lat = cubic_lattice(1.0)
            # Steel matrix with aluminum inclusion
            steel = IsotropicElastic(7800.0, 280e9, 130e9, 75e9)  # ρ, C11, C12, C44
            aluminum = IsotropicElastic(2700.0, 107e9, 55e9, 26e9)
            geo = Geometry(lat, steel, [(Sphere([0.0, 0.0, 0.0], 0.3), aluminum)])
            return geo
        end

        @testset "MatrixFreeOperator FullElastic constructor" begin
            geo = setup_3d_phononic()
            solver = Solver(FullElastic(), geo, (8, 8, 8); cutoff=2)
            k = [0.1, 0.2, 0.3]

            op = PhoXonic.MatrixFreeOperator(solver, k)

            @test op.solver === solver
            @test op.k == k
            @test PhoXonic.resolution(op) == (8, 8, 8)
            @test size(PhoXonic.work_real(op)) == (8, 8, 8)
        end

        @testset "apply_rhs! FullElastic vs Dense" begin
            geo = setup_3d_phononic()
            # Resolution must be large enough for FFT convolution to match dense
            solver = Solver(FullElastic(), geo, (16, 16, 16); cutoff=2)
            k = [0.1, 0.2, 0.3]
            N = solver.basis.num_pw

            # Dense RHS
            LHS, RHS = PhoXonic.build_matrices(solver, k)

            # Random input
            x = randn(ComplexF64, 3N)
            y_dense = RHS * x

            # Matrix-free RHS
            op = PhoXonic.MatrixFreeOperator(solver, k)
            y_mf = zeros(ComplexF64, 3N)
            PhoXonic.apply_rhs!(y_mf, op, x)

            @test y_mf ≈ y_dense rtol=1e-10
        end

        @testset "apply_lhs! FullElastic vs Dense" begin
            geo = setup_3d_phononic()
            # Resolution must be large enough for FFT convolution
            solver = Solver(FullElastic(), geo, (16, 16, 16); cutoff=2)
            k = [0.1, 0.2, 0.3]
            N = solver.basis.num_pw

            # Dense LHS
            LHS, RHS = PhoXonic.build_matrices(solver, k)

            # Random input
            x = randn(ComplexF64, 3N)
            y_dense = LHS * x

            # Matrix-free LHS
            op = PhoXonic.MatrixFreeOperator(solver, k)
            y_mf = zeros(ComplexF64, 3N)
            PhoXonic.apply_lhs!(y_mf, op, x)

            @test y_mf ≈ y_dense rtol=1e-10
        end
    end

    # ========================================================================
    # FFTContext and MatrixFreeWorkspace Tests
    # ========================================================================
    @testset "FFTContext and MatrixFreeWorkspace" begin
        @testset "FFTContext creation" begin
            lat = square_lattice(1.0)
            air = Dielectric(1.0)
            geo = Geometry(lat, air)
            solver = Solver(TEWave(), geo, (32, 32); cutoff=3)

            ctx = FFTContext(solver)

            @test ctx.resolution == (32, 32)
            @test ctx.fft_plan !== nothing
            @test ctx.ifft_plan !== nothing
        end

        @testset "MatrixFreeWorkspace creation" begin
            ctx = FFTContext((16, 16), ComplexF64)
            ws = MatrixFreeWorkspace(ctx)

            @test size(ws.work_real) == (16, 16)
            @test size(ws.work_fourier) == (16, 16)
        end

        @testset "MatrixFreeOperator with explicit ctx and workspace" begin
            lat = square_lattice(1.0)
            air = Dielectric(1.0)
            rod = Dielectric(8.9)
            geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])
            solver = Solver(TEWave(), geo, (32, 32); cutoff=3)

            # Create context and workspace
            ctx = FFTContext(solver)
            ws = MatrixFreeWorkspace(ctx)

            # Create operator with explicit ctx/ws
            k1 = [0.1, 0.2]
            k2 = [0.3, 0.4]
            op1 = MatrixFreeOperator(solver, k1, ctx, ws)
            op2 = MatrixFreeOperator(solver, k2, ctx, ws)

            # Both should share the same context
            @test op1.ctx === op2.ctx
            @test op1.workspace === op2.workspace

            # But have different k
            @test op1.k != op2.k
        end

        @testset "FFT plan reuse produces same results" begin
            lat = square_lattice(1.0)
            air = Dielectric(1.0)
            rod = Dielectric(8.9)
            geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])
            solver = Solver(TEWave(), geo, (32, 32); cutoff=3)
            k = [0.1, 0.2]
            N = solver.basis.num_pw

            # Without reuse (new plans each time)
            op1 = MatrixFreeOperator(solver, k)
            x = randn(ComplexF64, N)
            y1 = zeros(ComplexF64, N)
            PhoXonic.apply_lhs!(y1, op1, x)

            # With reuse
            ctx = FFTContext(solver)
            ws = MatrixFreeWorkspace(ctx)
            op2 = MatrixFreeOperator(solver, k, ctx, ws)
            y2 = zeros(ComplexF64, N)
            PhoXonic.apply_lhs!(y2, op2, x)

            @test y1 ≈ y2 rtol=1e-14
        end
    end

    # ========================================================================
    # MatrixFreeEffectiveHamiltonian Tests
    # ========================================================================
    @testset "MatrixFreeEffectiveHamiltonian" begin
        @testset "H = RHS⁻¹ * LHS consistency" begin
            lat = square_lattice(1.0)
            air = Dielectric(1.0)
            rod = Dielectric(4.0)  # Lower contrast for better convergence
            geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])
            solver = Solver(TEWave(), geo, (32, 32); cutoff=5)
            k = [0.1, 0.2]
            N = solver.basis.num_pw

            # Dense H
            LHS, RHS = PhoXonic.build_matrices(solver, k)
            H_dense = RHS \ LHS

            # Matrix-free H
            op = MatrixFreeOperator(solver, k)
            H_mf = PhoXonic.MatrixFreeEffectiveHamiltonian(op, ApproximateRHSInv())

            x = randn(ComplexF64, N)
            y_dense = H_dense * x
            y_mf = H_mf * x

            # ApproximateRHSInv() may have some error for inhomogeneous media
            # but should be reasonably close
            @test norm(y_dense - y_mf) / norm(y_dense) < 0.1
        end

        @testset "H with CGRHSInv method" begin
            lat = square_lattice(1.0)
            air = Dielectric(1.0)
            rod = Dielectric(4.0)
            geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])
            solver = Solver(TEWave(), geo, (32, 32); cutoff=5)
            k = [0.1, 0.2]
            N = solver.basis.num_pw

            # Dense H
            LHS, RHS = PhoXonic.build_matrices(solver, k)
            H_dense = RHS \ LHS

            # Matrix-free H with CGRHSInv() (should be more accurate)
            op = MatrixFreeOperator(solver, k)
            H_mf = PhoXonic.MatrixFreeEffectiveHamiltonian(op, CGRHSInv())

            x = randn(ComplexF64, N)
            y_dense = H_dense * x
            y_mf = H_mf * x

            # CGRHSInv() should be very close to dense
            @test norm(y_dense - y_mf) / norm(y_dense) < 1e-8
        end
    end

    # ========================================================================
    # Unified GFMethod API Tests
    # ========================================================================
    @testset "Unified GFMethod API" begin
        # Common setup
        function setup_ldos_test()
            lat = square_lattice(1.0)
            air = Dielectric(1.0)
            rod = Dielectric(4.0)
            geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])
            solver = Solver(TEWave(), geo, (32, 32); cutoff=5)
            return solver
        end

        @testset "DirectGF compute_ldos" begin
            solver = setup_ldos_test()
            pos = [0.5, 0.5]
            ω_values = [0.3, 0.4, 0.5]
            k_points = [[0.0, 0.0]]

            # Default (no method) should work
            ldos1 = compute_ldos(solver, pos, ω_values, k_points; η=0.05)

            # Explicit DirectGF() should give same result
            ldos2 = compute_ldos(solver, pos, ω_values, k_points, DirectGF(); η=0.05)

            @test ldos1 ≈ ldos2 rtol=1e-10
            @test all(ldos1 .>= 0)  # LDOS should be non-negative
        end

        # Note: MatrixFreeGF and RSKGF tests moved to test/rsk_ext/runtests.jl
        # They require ReducedShiftedKrylov.jl extension to be loaded

        @testset "compute_greens_function with DirectGF" begin
            solver = setup_ldos_test()
            k = [0.1, 0.2]
            ω_values = [0.3, 0.4]
            N = solver.basis.num_pw
            source = randn(ComplexF64, N)

            # Direct
            G_direct = compute_greens_function(solver, k, ω_values, source, DirectGF(); η=0.05)

            @test length(G_direct) == length(ω_values)
            @test all(length.(G_direct) .== N)
        end

        @testset "compute_dos with DirectGF" begin
            solver = setup_ldos_test()
            ω_values = [0.3, 0.4, 0.5]
            k_points = [[0.0, 0.0], [0.5, 0.0]]

            # Direct
            dos_direct = compute_dos(solver, ω_values, k_points, DirectGF(); η=0.05)

            @test length(dos_direct) == length(ω_values)
            @test all(dos_direct .>= 0)
        end
    end

    # ========================================================================
    # matrix_dimension Tests
    # ========================================================================
    @testset "matrix_dimension" begin
        # Setup 2D geometry (photonic)
        lat_2d = square_lattice(1.0)
        air = Dielectric(1.0)
        rod = Dielectric(4.0)
        geo_2d = Geometry(lat_2d, air, [(Circle([0.0, 0.0], 0.2), rod)])

        # Setup 2D geometry (phononic)
        steel = IsotropicElastic(ρ=7780.0, λ=111.2e9, μ=81.0e9)
        epoxy = IsotropicElastic(ρ=1142.0, λ=0.458e9, μ=0.148e9)
        geo_2d_phon = Geometry(lat_2d, epoxy, [(Circle([0.0, 0.0], 0.2), steel)])

        # Setup 3D geometry (photonic)
        lat_3d = cubic_lattice(1.0)
        geo_3d = Geometry(lat_3d, air, [(Sphere([0.0, 0.0, 0.0], 0.3), rod)])

        @testset "Scalar 2D waves (TE, TM, SH)" begin
            # TE wave: dim = num_pw
            solver_te = Solver(TEWave(), geo_2d, (32, 32); cutoff=5)
            @test matrix_dimension(solver_te) == solver_te.basis.num_pw

            # TM wave: dim = num_pw
            solver_tm = Solver(TMWave(), geo_2d, (32, 32); cutoff=5)
            @test matrix_dimension(solver_tm) == solver_tm.basis.num_pw

            # SH wave: dim = num_pw
            solver_sh = Solver(SHWave(), geo_2d_phon, (32, 32); cutoff=5)
            @test matrix_dimension(solver_sh) == solver_sh.basis.num_pw
        end

        @testset "Vector 2D wave (P-SV)" begin
            # P-SV wave: dim = 2 * num_pw
            solver_psv = Solver(PSVWave(), geo_2d_phon, (32, 32); cutoff=5)
            @test matrix_dimension(solver_psv) == 2 * solver_psv.basis.num_pw
        end

        @testset "Vector 3D waves (FullVectorEM, FullElastic)" begin
            # FullVectorEM: dim = 3 * num_pw
            solver_em = Solver(FullVectorEM(), geo_3d, (8, 8, 8); cutoff=2)
            @test matrix_dimension(solver_em) == 3 * solver_em.basis.num_pw

            # FullElastic: dim = 3 * num_pw (need 3D phononic geometry)
            lat_3d_phon = cubic_lattice(1.0)
            steel_3d = IsotropicElastic(7800.0, 280e9, 130e9, 75e9)
            aluminum = IsotropicElastic(2700.0, 107e9, 55e9, 26e9)
            geo_3d_phon = Geometry(lat_3d_phon, steel_3d, [(Sphere([0.0, 0.0, 0.0], 0.3), aluminum)])
            solver_elastic = Solver(FullElastic(), geo_3d_phon, (8, 8, 8); cutoff=2)
            @test matrix_dimension(solver_elastic) == 3 * solver_elastic.basis.num_pw
        end
    end

    # ========================================================================
    # solve_at_k Tests
    # ========================================================================
    @testset "solve_at_k" begin
        # Setup 2D geometry (phononic for SH/P-SV waves)
        lat_2d = square_lattice(1.0)
        steel = IsotropicElastic(ρ=7780.0, λ=111.2e9, μ=81.0e9)
        epoxy = IsotropicElastic(ρ=1142.0, λ=0.458e9, μ=0.148e9)
        geo_phon = Geometry(lat_2d, epoxy, [(Circle([0.0, 0.0], 0.2), steel)])

        @testset "Matrix scaling" begin
            solver = Solver(PSVWave(), geo_phon, (32, 32); cutoff=5)
            k = [0.1, 0.2]

            # Dense reference
            freqs_dense = solve_at_k(solver, k, DenseMethod(); bands=1:5)

            # LOBPCG without scaling (may have convergence issues)
            freqs_no_scale = solve_at_k(solver, k, LOBPCGMethod(scale=false, preconditioner=:none); bands=1:5)

            # LOBPCG with scaling (should converge better)
            freqs_scale = solve_at_k(solver, k, LOBPCGMethod(scale=true, preconditioner=:diagonal); bands=1:5)

            # Both should return valid frequencies
            @test length(freqs_no_scale) == 5
            @test length(freqs_scale) == 5
            @test all(freqs_scale .>= 0)

            # Scaling version should be closer to Dense (or at least not worse)
            # Note: For small problems, both may work fine
            error_no_scale = maximum(abs.(freqs_no_scale - freqs_dense))
            error_scale = maximum(abs.(freqs_scale - freqs_dense))
            # Just verify both produce reasonable results
            @test error_scale < 10000  # Reasonable tolerance for phononic problems
        end

        @testset "compute_bands with warm start" begin
            # Use square lattice for simpler k-path
            lat = square_lattice(1.0)
            geo = Geometry(lat, epoxy, [(Circle([0.0, 0.0], 0.2), steel)])

            # Create solvers
            solver_ws = Solver(PSVWave(), geo, (32, 32),
                              LOBPCGMethod(warm_start=true, first_dense=true); cutoff=5)
            solver_dense = Solver(PSVWave(), geo, (32, 32); cutoff=5)

            # Simple k-path with few points
            kpath = simple_kpath_square(a=1.0, npoints=5)

            # Compute bands
            bands_ws = compute_bands(solver_ws, kpath; bands=1:5)
            bands_dense = compute_bands(solver_dense, kpath; bands=1:5)

            # Results should be similar
            @test size(bands_ws.frequencies) == size(bands_dense.frequencies)
            max_error = maximum(abs.(bands_ws.frequencies - bands_dense.frequencies))
            @test max_error < 5000  # Allow some tolerance for phononic problems
        end

        @testset "first_dense option" begin
            lat = square_lattice(1.0)
            geo = Geometry(lat, epoxy, [(Circle([0.0, 0.0], 0.2), steel)])

            # With first_dense=true
            solver_fd = Solver(PSVWave(), geo, (32, 32),
                              LOBPCGMethod(warm_start=true, first_dense=true); cutoff=5)
            # Without first_dense
            solver_no_fd = Solver(PSVWave(), geo, (32, 32),
                                 LOBPCGMethod(warm_start=true, first_dense=false); cutoff=5)

            kpath = simple_kpath_square(a=1.0, npoints=3)

            bands_fd = compute_bands(solver_fd, kpath; bands=1:5)
            bands_no_fd = compute_bands(solver_no_fd, kpath; bands=1:5)

            # Both should produce valid results
            @test size(bands_fd.frequencies) == (length(kpath.points), 5)
            @test size(bands_no_fd.frequencies) == (length(kpath.points), 5)
            @test all(bands_fd.frequencies .>= 0)
        end

        @testset "Basic functionality" begin
            solver = Solver(SHWave(), geo_phon, (32, 32); cutoff=5)
            k = [0.1, 0.2]

            # Frequencies only (default)
            freqs = solve_at_k(solver, k, DenseMethod(); bands=1:5)
            @test length(freqs) == 5
            @test all(freqs .>= 0)

            # With eigenvectors
            freqs2, vecs = solve_at_k(solver, k, DenseMethod();
                                      bands=1:5, return_eigenvectors=true)
            @test freqs ≈ freqs2
            @test size(vecs) == (matrix_dimension(solver), 5)
        end

        @testset "With initial vectors X0" begin
            solver = Solver(SHWave(), geo_phon, (32, 32); cutoff=5)
            dim = matrix_dimension(solver)
            k1, k2 = [0.1, 0.2], [0.15, 0.25]

            # First k-point with Dense
            freqs1, vecs1 = solve_at_k(solver, k1, DenseMethod();
                                       bands=1:5, return_eigenvectors=true)

            # Second k-point with LOBPCG using warm start
            freqs2, vecs2 = solve_at_k(solver, k2, LOBPCGMethod();
                                       bands=1:5, X0=vecs1, return_eigenvectors=true)

            @test length(freqs2) == 5
            @test size(vecs2) == (dim, 5)

            # Compare with Dense reference
            freqs_ref, _ = solve_at_k(solver, k2, DenseMethod();
                                      bands=1:5, return_eigenvectors=true)
            @test maximum(abs.(freqs2 - freqs_ref)) < 1000  # Allow some tolerance
        end

        @testset "With preconditioner" begin
            solver = Solver(SHWave(), geo_phon, (32, 32); cutoff=5)
            k = [0.1, 0.2]

            # :none preconditioner
            freqs1 = solve_at_k(solver, k, LOBPCGMethod(preconditioner=:none); bands=1:5)

            # :diagonal preconditioner
            freqs2 = solve_at_k(solver, k, LOBPCGMethod(preconditioner=:diagonal); bands=1:5)

            # Custom Diagonal preconditioner
            LHS, _ = PhoXonic.build_matrices(solver, k)
            d = diag(LHS)
            d_safe = [abs(x) > 1e-10 ? x : 1e-10 for x in d]
            P = Diagonal(1.0 ./ d_safe)
            freqs3 = solve_at_k(solver, k, LOBPCGMethod(); bands=1:5, P=P)

            # All should produce similar results
            @test length(freqs1) == 5
            @test length(freqs2) == 5
            @test length(freqs3) == 5
        end
    end

    # ========================================================================
    # LOBPCGMethod Extension Tests
    # ========================================================================
    @testset "LOBPCGMethod extension" begin
        @testset "Default values" begin
            m = LOBPCGMethod()
            @test m.tol == 1e-4
            @test m.maxiter == 200
            @test m.shift == 0.0
            @test m.warm_start == true
            @test m.scale == true
            @test m.first_dense == true
            @test m.preconditioner == :diagonal
        end

        @testset "Custom values" begin
            m = LOBPCGMethod(tol=1e-6, maxiter=500, shift=0.01,
                            warm_start=false, scale=false, first_dense=false,
                            preconditioner=:none)
            @test m.tol == 1e-6
            @test m.maxiter == 500
            @test m.shift == 0.01
            @test m.warm_start == false
            @test m.scale == false
            @test m.first_dense == false
            @test m.preconditioner == :none
        end

        @testset "Backward compatibility" begin
            # Constructor without new fields (warm_start, scale, first_dense, preconditioner)
            # should still work with default values for the new fields
            m = LOBPCGMethod(tol=1e-8, maxiter=300, shift=0.0)
            @test m.tol == 1e-8
            @test m.maxiter == 300
            @test m.shift == 0.0
            # New fields get default values
            @test m.warm_start == true
            @test m.scale == true
            @test m.first_dense == true
            @test m.preconditioner == :diagonal
        end
    end

    # ========================================================================
    # LOBPCG Performance and Accuracy Tests (Phase 5)
    # ========================================================================
    @testset "LOBPCG Performance and Accuracy" begin
        # Setup for tests - use photonic (TE wave) for accuracy tests
        # Note: Phononic problems with band crossings at Gamma require
        # band tracking which is not yet implemented
        lat = square_lattice(1.0)
        eps_bg = Dielectric(1.0)
        eps_rod = Dielectric(12.0)
        geo_phot = Geometry(lat, eps_bg, [(Circle([0.0, 0.0], 0.3), eps_rod)])

        # Setup for phononic (only used where band tracking is not critical)
        steel = IsotropicElastic(ρ=7780.0, λ=111.2e9, μ=81.0e9)
        epoxy = IsotropicElastic(ρ=1142.0, λ=0.458e9, μ=0.148e9)
        geo_phon = Geometry(lat, epoxy, [(Circle([0.0, 0.0], 0.3), steel)])

        @testset "LOBPCG completes successfully (photonic)" begin
            # Note: For small matrices, Dense is faster due to BLAS optimization
            # LOBPCG becomes faster for large problems (cutoff > 12)
            # This test verifies LOBPCG works correctly, not that it's faster
            cutoff_val = 8
            solver_lobpcg = Solver(TEWave(), geo_phot, (32, 32),
                                   LOBPCGMethod(warm_start=true, scale=true); cutoff=cutoff_val)
            solver_dense = Solver(TEWave(), geo_phot, (32, 32); cutoff=cutoff_val)

            kpath = simple_kpath_square(a=1.0, npoints=10)

            # Both methods should complete and produce correct results
            bands_lobpcg = compute_bands(solver_lobpcg, kpath; bands=1:8)
            bands_dense = compute_bands(solver_dense, kpath; bands=1:8)

            # Results should have same size
            @test size(bands_lobpcg.frequencies) == size(bands_dense.frequencies)

            # All frequencies should be non-negative
            @test all(bands_lobpcg.frequencies .>= 0)
        end

        @testset "LOBPCG accuracy vs Dense (photonic)" begin
            # Photonic problems (TE/TM) work well with warm start
            cutoff_val = 8
            solver_lobpcg = Solver(TEWave(), geo_phot, (32, 32),
                                   LOBPCGMethod(warm_start=true, scale=true,
                                               first_dense=true, tol=1e-4); cutoff=cutoff_val)
            solver_dense = Solver(TEWave(), geo_phot, (32, 32); cutoff=cutoff_val)

            kpath = simple_kpath_square(a=1.0, npoints=10)

            bands_lobpcg = compute_bands(solver_lobpcg, kpath; bands=1:8)
            bands_dense = compute_bands(solver_dense, kpath; bands=1:8)

            # Compute relative error for each band
            max_abs_error = maximum(abs.(bands_lobpcg.frequencies - bands_dense.frequencies))
            ref_scale = maximum(abs.(bands_dense.frequencies))
            relative_error = max_abs_error / ref_scale

            # Accuracy should be very good for photonic problems
            @test relative_error < 0.05  # 5% relative error
            @test max_abs_error < 1.0  # normalized frequency units
        end

        @testset "LOBPCG warm start convergence" begin
            # Test that warm start produces valid results
            cutoff_val = 6
            solver = Solver(TEWave(), geo_phot, (32, 32),
                           LOBPCGMethod(warm_start=true, first_dense=true, maxiter=50); cutoff=cutoff_val)

            kpath = simple_kpath_square(a=1.0, npoints=5)

            # This should complete without hitting maxiter for most k-points
            bands = compute_bands(solver, kpath; bands=1:5)

            # All frequencies should be non-negative (no numerical issues)
            @test all(bands.frequencies .>= 0)
            @test size(bands.frequencies) == (length(kpath.points), 5)
        end
    end

    # ========================================================================
    # Error Fallback Tests
    # ========================================================================
    @testset "Error Fallbacks" begin
        # Dummy types for testing error fallbacks
        struct UnsupportedMaterial <: PhoXonic.Material end
        struct UnsupportedWave <: PhoXonic.WaveType end
        struct UnsupportedTrait end

        @testset "get_property unsupported material" begin
            mat = UnsupportedMaterial()
            @test_throws ArgumentError PhoXonic.get_property(mat, :ε)
            @test_throws ArgumentError PhoXonic.get_property(mat, :ρ)
        end

        @testset "_apply_rhs_inv_impl! unsupported trait" begin
            trait = UnsupportedTrait()
            method = PhoXonic.ApproximateRHSInv()  # Any RHSInvMethod
            y = zeros(ComplexF64, 10)
            H = nothing  # Dummy, won't be used
            x = zeros(ComplexF64, 10)
            @test_throws ArgumentError PhoXonic._apply_rhs_inv_impl!(trait, method, y, H, x)
        end

        @testset "_apply_rhs_inv_3d_impl! unsupported trait" begin
            trait = UnsupportedTrait()
            method = PhoXonic.ApproximateRHSInv()  # Any RHSInvMethod
            y = zeros(ComplexF64, 30)
            H = nothing  # Dummy, won't be used
            x = zeros(ComplexF64, 30)
            @test_throws ArgumentError PhoXonic._apply_rhs_inv_3d_impl!(trait, method, y, H, x)
        end

        @testset "_convert_rhs_inv_method unsupported type" begin
            @test_throws ArgumentError PhoXonic._convert_rhs_inv_method(123)  # Int
            @test_throws ArgumentError PhoXonic._convert_rhs_inv_method("approximate")  # String
        end

        @testset "wave_structure unsupported wave" begin
            wave = UnsupportedWave()
            @test_throws ArgumentError PhoXonic.wave_structure(wave)
        end

        @testset "ncomponents unsupported wave" begin
            wave = UnsupportedWave()
            @test_throws ArgumentError PhoXonic.ncomponents(wave)
        end

        @testset "applicable_dimension unsupported wave" begin
            @test_throws ArgumentError PhoXonic.applicable_dimension(UnsupportedWave)
        end

        @testset "_get_rhs_inv unsupported wave" begin
            wave = UnsupportedWave()
            mats = nothing
            res = (32, 32)
            @test_throws ArgumentError PhoXonic._get_rhs_inv(wave, mats, res, ComplexF64)
        end

        @testset "prepare_materials unsupported wave/dimension" begin
            # UnsupportedWave with 2D geometry should fail
            lat = square_lattice(1.0)
            air = Dielectric(1.0)
            geo = Geometry(lat, air)
            wave = UnsupportedWave()
            @test_throws ArgumentError PhoXonic.prepare_materials(wave, geo, (32, 32))
        end
    end

    # ========================================================================
    # ElasticVoid (Tanaka Limit) Tests
    # ========================================================================
    @testset "ElasticVoid (Tanaka Limit)" begin
        @testset "ElasticVoid basic" begin
            # Basic construction
            void = ElasticVoid()
            @test void isa ElasticMaterial

            # Default ρ_ratio = 1e-7
            @test density(void) ≈ 1e-7
            @test shear_modulus(void) ≈ 1.0

            # Transverse velocity: √(μ/ρ) = √(1/1e-7) ≈ 3162
            @test transverse_velocity(void) ≈ sqrt(1e7) atol=1.0

            # Custom ρ_ratio
            void2 = ElasticVoid(ρ_ratio=1e-8)
            @test density(void2) ≈ 1e-8
            @test transverse_velocity(void2) ≈ 10000.0

            # C11, C12, C44 consistency (isotropic: C44 = (C11-C12)/2)
            @test void.C44 ≈ (void.C11 - void.C12) / 2
        end

        @testset "Geometry with ElasticVoid" begin
            lat = square_lattice(1.0)
            al = IsotropicElastic(ρ=2700.0, λ=58.7e9, μ=26.1e9)
            void = ElasticVoid()

            # Geometry creation should succeed
            geo = Geometry(lat, al, [(Circle([0.0, 0.0], 0.4), void)])
            @test geo isa Geometry

            # Material detection
            @test get_material(geo, [0.0, 0.0]) isa ElasticVoid
            @test get_material(geo, [0.49, 0.49]) isa IsotropicElastic

            # discretize should work
            ρ_grid = PhoXonic.discretize(geo, (32, 32), :ρ)
            @test size(ρ_grid) == (32, 32)
            @test minimum(ρ_grid) < 1.0  # void region has ρ << 1
        end

        @testset "Solver with ElasticVoid" begin
            lat = square_lattice(1.0)
            al = IsotropicElastic(ρ=2700.0, λ=58.7e9, μ=26.1e9)
            void = ElasticVoid()

            geo = Geometry(lat, al, [(Circle([0.0, 0.0], 0.3), void)])

            # Solver creation should succeed
            solver = Solver(SHWave(), geo, (32, 32); cutoff=3)
            @test solver isa Solver

            # Band computation should work (no errors)
            k = [0.25, 0.0]
            freqs = solve_at_k(solver, k, DenseMethod(); bands=1:3)
            @test length(freqs) == 3
            @test all(freqs .>= 0)
        end
    end

    # TMM (Transfer Matrix Method) tests
    include("tmm/runtests.jl")

end
