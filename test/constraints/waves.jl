@testitem "Wave constructors" begin
    @test AboveSin([0, 0, 1], [1, 0, 0], 5.0, 2.0, 4.0, 0.0) ==
          Wave{Over,Float64}([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 5.0, 2.0, 4.0, 0.0, 5.0)
    @test BelowSin(up=[0, 0, 1], along=[1, 0, 0], d0=5.0, amplitude=2.0, wavelength=4.0, phase=0.0) ==
          Wave{Below,Float64}([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 5.0, 2.0, 4.0, 0.0, 5.0)
    # cos is sin shifted by pi/2
    @test AboveCos([0, 0, 1], [1, 0, 0], 5.0, 2.0, 4.0, 0.0) ==
          Wave{Over,Float64}([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 5.0, 2.0, 4.0, pi / 2, 5.0)
    @test BelowCos(up=[0, 0, 1], along=[1, 0, 0], d0=5.0, amplitude=2.0, wavelength=4.0, phase=0.0) ==
          Wave{Below,Float64}([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 5.0, 2.0, 4.0, pi / 2, 5.0)

    @test AboveRadialSin([0, 0, 1], [0, 0, 0], 5.0, 2.0, 4.0, 0.0) ==
          RadialWave{Over,Float64}([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 5.0, 2.0, 4.0, 0.0, 5.0)
    @test AboveRadialCos([0, 0, 1], [0, 0, 0], 5.0, 2.0, 4.0, 0.0) ==
          RadialWave{Over,Float64}([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 5.0, 2.0, 4.0, pi / 2, 5.0)
end

@testitem "Wave degenerates to Plane at amplitude=0" begin
    using StaticArrays
    x = SVector{3,Float64}(1.0, 2.0, 3.0)
    w = AboveSin([0, 0, 1], [1, 0, 0], 5.0, 0.0, 4.0, 0.3)
    p = AbovePlane([0, 0, 1], 5.0)
    @test Packmol.constraint_penalty(w, x) ≈ Packmol.constraint_penalty(p, x)
    @test Packmol.constraint_gradient(w, x) ≈ Packmol.constraint_gradient(p, x)
end

@testitem "Wave gradients" begin
    using ForwardDiff
    using StaticArrays
    points = [
        SVector{3,Float64}(1.0, 2.0, 3.0),
        SVector{3,Float64}(0.0, 0.0, 10.0),
        SVector{3,Float64}(-2.0, 3.0, -1.0),
        SVector{3,Float64}(0.5, -1.5, 0.2),
    ]
    for c in (
        AboveSin([0, 0, 1], [1, 0, 0], 5.0, 2.0, 4.0, 0.3),
        BelowSin([0, 0, 1], [1, 0, 0], 5.0, 2.0, 4.0, 0.3),
        AboveCos([1, 1, 1] ./ sqrt(3.0), [1, -1, 0] ./ sqrt(2.0), 1.0, -3.0, 2.5, 0.1),
        BelowCos([1, 1, 1] ./ sqrt(3.0), [1, -1, 0] ./ sqrt(2.0), 1.0, -3.0, 2.5, 0.1),
    )
        for x in points
            @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
        end
    end
end

@testitem "Wave input file parsing" begin
    structure_data = Dict{Symbol,Any}(:filename => "dummy.pdb")
    data = split("0. 0. 1.  1. 0. 0.  5.0 2.0 4.0 0.3")
    @test Packmol.parse_constraint["above sin"](structure_data, data) ==
          AboveSin([0, 0, 1], [1, 0, 0], 5.0, 2.0, 4.0, 0.3)
    @test Packmol.parse_constraint["below sin"](structure_data, data) ==
          BelowSin([0, 0, 1], [1, 0, 0], 5.0, 2.0, 4.0, 0.3)
    @test Packmol.parse_constraint["above cos"](structure_data, data) ==
          AboveCos([0, 0, 1], [1, 0, 0], 5.0, 2.0, 4.0, 0.3)
    @test Packmol.parse_constraint["below cos"](structure_data, data) ==
          BelowCos([0, 0, 1], [1, 0, 0], 5.0, 2.0, 4.0, 0.3)

    data = split("0. 0. 1.  0. 0. 0.  5.0 2.0 4.0 0.3")
    @test Packmol.parse_constraint["above radial_sin"](structure_data, data) ==
          AboveRadialSin([0, 0, 1], [0, 0, 0], 5.0, 2.0, 4.0, 0.3)
    @test Packmol.parse_constraint["above radial_cos"](structure_data, data) ==
          AboveRadialCos([0, 0, 1], [0, 0, 0], 5.0, 2.0, 4.0, 0.3)
end

@testitem "RadialWave gradients" begin
    using ForwardDiff
    using StaticArrays
    # Away from r = 0, the surface is smooth: compare against ForwardDiff.
    points = [
        SVector{3,Float64}(1.0, 2.0, 3.0),
        SVector{3,Float64}(0.0, 0.0, 10.0),
        SVector{3,Float64}(-2.0, 3.0, -1.0),
        SVector{3,Float64}(1e-3, 0.0, 0.0),
    ]
    for c in (
        AboveRadialSin([0, 0, 1], [0, 0, 0], 5.0, 2.0, 4.0, 0.3),
        BelowRadialSin([0, 0, 1], [0, 0, 0], 5.0, 2.0, 4.0, 0.3),
        AboveRadialCos([1, 1, 1] ./ sqrt(3.0), [0.2, -0.3, 0.1], 1.0, -3.0, 2.5, 0.1),
    )
        for x in points
            @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
        end
    end
    # At r = 0 (the ring's center), the gradient's radial contribution is a
    # removable singularity handled by treating it as zero: check the
    # analytic gradient is finite there (not NaN/Inf from a 0/0 division).
    c = AboveRadialSin([0, 0, 1], [0, 0, 0], 5.0, 2.0, 4.0, 0.3)
    x = SVector{3,Float64}(0.0, 0.0, 0.0)
    g = Packmol.constraint_gradient(c, x)
    @test all(isfinite, g)
end
