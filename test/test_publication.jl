using Test
using PhoXonic
using PythonPlot
using StaticArrays

@testset "savefig_publication (PythonPlot)" begin
    # Lightweight 2D BandStructure fixture (4 k-points × 3 bands)
    bs2 = PhoXonic.BandStructure{2}(
        [SVector(0.0, 0.0), SVector(0.5, 0.0), SVector(0.5, 0.5), SVector(0.0, 0.0)],
        [0.0, 1.0, 2.0, 3.0],
        Float64[0.10 0.30 0.55;
                0.20 0.40 0.65;
                0.15 0.35 0.60;
                0.05 0.25 0.50],
        [(1, "Γ"), (2, "X"), (3, "M"), (4, "Γ")],
    )

    # Lightweight 3D BandStructure fixture
    bs3 = PhoXonic.BandStructure{3}(
        [SVector(0.0, 0.0, 0.0), SVector(0.5, 0.0, 0.0),
         SVector(0.5, 0.5, 0.0), SVector(0.5, 0.5, 0.5),
         SVector(0.0, 0.0, 0.0)],
        [0.0, 1.0, 2.0, 3.0, 4.0],
        Float64[0.10 0.30 0.55;
                0.20 0.40 0.65;
                0.15 0.35 0.60;
                0.18 0.38 0.62;
                0.05 0.25 0.50],
        [(1, "Γ"), (2, "X"), (3, "M"), (4, "R"), (5, "Γ")],
    )

    @testset "BandStructure{2} 単体" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bs2.pdf")
            ret = savefig_publication(bs2, path)
            @test isfile(path)
            @test ret == path
            @test filesize(path) > 0
        end
    end

    @testset "BandStructure{3} 単体" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bs3.pdf")
            savefig_publication(bs3, path)
            @test isfile(path)
        end
    end

    @testset "axis_width_cm / axis_height_cm カスタム" begin
        mktempdir() do tmp
            path = joinpath(tmp, "custom_size.pdf")
            savefig_publication(bs2, path; axis_width_cm=10.0, axis_height_cm=5.0)
            @test isfile(path)
        end
    end

    @testset "title kwarg" begin
        mktempdir() do tmp
            path = joinpath(tmp, "with_title.pdf")
            savefig_publication(bs2, path; title="Test Bands")
            @test isfile(path)
        end
    end

    @testset "ylims kwarg" begin
        mktempdir() do tmp
            path = joinpath(tmp, "ylims.pdf")
            savefig_publication(bs2, path; ylims=(0.0, 1.0))
            @test isfile(path)
        end
    end

    @testset "show_gaps" begin
        mktempdir() do tmp
            path = joinpath(tmp, "gaps.pdf")
            savefig_publication(bs2, path; show_gaps=true)
            @test isfile(path)
        end
    end

    @testset "subplot layout=(1,2)" begin
        mktempdir() do tmp
            path = joinpath(tmp, "subplot_1x2.pdf")
            savefig_publication([bs2, bs2], path; layout=(1, 2))
            @test isfile(path)
        end
    end

    @testset "subplot layout=(2,1)" begin
        mktempdir() do tmp
            path = joinpath(tmp, "subplot_2x1.pdf")
            savefig_publication([bs2, bs2], path; layout=(2, 1))
            @test isfile(path)
        end
    end

    @testset "overlay モード" begin
        mktempdir() do tmp
            path = joinpath(tmp, "overlay.pdf")
            savefig_publication([bs2, bs2], path; overlay=true)
            @test isfile(path)
        end
    end

    @testset "PNG 出力" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bs2.png")
            savefig_publication(bs2, path)
            @test isfile(path)
        end
    end
end
