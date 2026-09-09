@testitem "Gaussian constructors" begin
    @test AboveGaussian([0, 0, 1], [1, 0, 0], 5.0, 2.0, 0.0, 1.0) ==
          Gaussian{Over,Float64}([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 5.0, 2.0, 0.0, 1.0, 5.0)
    @test AboveGaussian(up=[0, 0, 1], along=[1, 0, 0], d0=5.0, amplitude=2.0, center=0.0, sigma=1.0) ==
          Gaussian{Over,Float64}([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 5.0, 2.0, 0.0, 1.0, 5.0)
    @test BelowGaussian([0, 0, 1], [1, 0, 0], 5.0, 2.0, 0.0, 1.0, 2.0) ==
          Gaussian{Below,Float64}([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 5.0, 2.0, 0.0, 1.0, 2.0)

    @test AboveRadialGaussian([0, 0, 1], [0, 0, 0], 5.0, 2.0, 1.0) ==
          RadialGaussian{Over,Float64}([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 5.0, 2.0, 1.0, 5.0)
    @test BelowRadialGaussian(up=[0, 0, 1], center=[0, 0, 0], d0=5.0, amplitude=2.0, sigma=1.0) ==
          RadialGaussian{Below,Float64}([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 5.0, 2.0, 1.0, 5.0)
end

@testitem "Gaussian degenerates to Plane at amplitude=0" begin
    using StaticArrays
    x = SVector{3,Float64}(1.0, 2.0, 3.0)
    g = AboveGaussian([0, 0, 1], [1, 0, 0], 5.0, 0.0, 0.0, 1.0)
    p = AbovePlane([0, 0, 1], 5.0)
    @test Packmol.constraint_penalty(g, x) ≈ Packmol.constraint_penalty(p, x)
    @test Packmol.constraint_gradient(g, x) ≈ Packmol.constraint_gradient(p, x)
end

@testitem "Gaussian gradients" begin
    using ForwardDiff
    using StaticArrays
    points = [
        SVector{3,Float64}(1.0, 2.0, 3.0),
        SVector{3,Float64}(0.0, 0.0, 10.0),
        SVector{3,Float64}(-2.0, 3.0, -1.0),
        SVector{3,Float64}(0.5, -1.5, 0.2),
    ]
    for c in (
        AboveGaussian([0, 0, 1], [1, 0, 0], 5.0, 2.0, 0.0, 1.5),
        BelowGaussian([0, 0, 1], [1, 0, 0], 5.0, 2.0, 0.0, 1.5),
        AboveGaussian([1, 1, 1] ./ sqrt(3.0), [1, -1, 0] ./ sqrt(2.0), 1.0, -3.0, 0.5, 0.7),
    )
        for x in points
            @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
        end
    end
end

@testitem "Gaussian input file parsing" begin
    structure_data = Dict{Symbol,Any}(:filename => "dummy.pdb")
    data = split("0. 0. 1.  1. 0. 0.  5.0 2.0 0.0 1.0")
    @test Packmol.parse_constraint["above gaussian"](structure_data, data) ==
          AboveGaussian([0, 0, 1], [1, 0, 0], 5.0, 2.0, 0.0, 1.0)
    @test Packmol.parse_constraint["below gaussian"](structure_data, data) ==
          BelowGaussian([0, 0, 1], [1, 0, 0], 5.0, 2.0, 0.0, 1.0)

    data = split("0. 0. 1.  0. 0. 0.  5.0 2.0 1.0")
    @test Packmol.parse_constraint["above radial_gaussian"](structure_data, data) ==
          AboveRadialGaussian([0, 0, 1], [0, 0, 0], 5.0, 2.0, 1.0)
    @test Packmol.parse_constraint["below radial_gaussian"](structure_data, data) ==
          BelowRadialGaussian([0, 0, 1], [0, 0, 0], 5.0, 2.0, 1.0)
end

@testitem "RadialGaussian gradients" begin
    using ForwardDiff
    using StaticArrays
    points = [
        SVector{3,Float64}(1.0, 2.0, 3.0),
        SVector{3,Float64}(0.0, 0.0, 10.0),
        SVector{3,Float64}(-2.0, 3.0, -1.0),
        SVector{3,Float64}(0.0, 0.0, 0.0), # exactly on the axis (r = 0)
    ]
    for c in (
        AboveRadialGaussian([0, 0, 1], [0, 0, 0], 5.0, 2.0, 1.5),
        BelowRadialGaussian([0, 0, 1], [0, 0, 0], 5.0, 2.0, 1.5),
        AboveRadialGaussian([1, 1, 1] ./ sqrt(3.0), [0.2, -0.3, 0.1], 1.0, -3.0, 0.7),
    )
        for x in points
            @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
        end
    end
end
