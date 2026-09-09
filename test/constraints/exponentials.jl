@testitem "Exponential constructors" begin
    @test AboveExponential([0, 0, 1], [1, 0, 0], 5.0, 2.0, 0.0, 0.5) ==
          Exponential{Over,Float64}([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 5.0, 2.0, 0.0, 0.5, 5.0)
    @test BelowExponential(up=[0, 0, 1], along=[1, 0, 0], d0=5.0, amplitude=2.0, center=0.0, rate=0.5) ==
          Exponential{Below,Float64}([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], 5.0, 2.0, 0.0, 0.5, 5.0)

    @test AboveRadialExponential([0, 0, 1], [0, 0, 0], 5.0, 2.0, 0.5) ==
          RadialExponential{Over,Float64}([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 5.0, 2.0, 0.5, 5.0)
    @test BelowRadialExponential(up=[0, 0, 1], center=[0, 0, 0], d0=5.0, amplitude=2.0, rate=0.5) ==
          RadialExponential{Below,Float64}([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 5.0, 2.0, 0.5, 5.0)
end

@testitem "Exponential degenerates to Plane at amplitude=0" begin
    using StaticArrays
    x = SVector{3,Float64}(1.0, 2.0, 3.0)
    e = AboveExponential([0, 0, 1], [1, 0, 0], 5.0, 0.0, 0.0, 0.5)
    p = AbovePlane([0, 0, 1], 5.0)
    @test Packmol.constraint_penalty(e, x) ≈ Packmol.constraint_penalty(p, x)
    @test Packmol.constraint_gradient(e, x) ≈ Packmol.constraint_gradient(p, x)
end

@testitem "Exponential gradients" begin
    using ForwardDiff
    using StaticArrays
    points = [
        SVector{3,Float64}(1.0, 2.0, 3.0),
        SVector{3,Float64}(0.0, 0.0, 10.0),
        SVector{3,Float64}(-2.0, 3.0, -1.0),
        SVector{3,Float64}(0.5, -1.5, 0.2),
    ]
    for c in (
        AboveExponential([0, 0, 1], [1, 0, 0], 5.0, 1.0, 0.0, 0.5),
        BelowExponential([0, 0, 1], [1, 0, 0], 5.0, 1.0, 0.0, 0.5),
        AboveExponential([1, 1, 1] ./ sqrt(3.0), [1, -1, 0] ./ sqrt(2.0), 1.0, -2.0, 0.3, -0.8),
    )
        for x in points
            @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
        end
    end
end

@testitem "Exponential input file parsing" begin
    structure_data = Dict{Symbol,Any}(:filename => "dummy.pdb")
    data = split("0. 0. 1.  1. 0. 0.  5.0 2.0 0.0 0.5")
    @test Packmol.parse_constraint["above exponential"](structure_data, data) ==
          AboveExponential([0, 0, 1], [1, 0, 0], 5.0, 2.0, 0.0, 0.5)
    @test Packmol.parse_constraint["below exponential"](structure_data, data) ==
          BelowExponential([0, 0, 1], [1, 0, 0], 5.0, 2.0, 0.0, 0.5)

    data = split("0. 0. 1.  0. 0. 0.  5.0 2.0 0.5")
    @test Packmol.parse_constraint["above radial_exponential"](structure_data, data) ==
          AboveRadialExponential([0, 0, 1], [0, 0, 0], 5.0, 2.0, 0.5)
    @test Packmol.parse_constraint["below radial_exponential"](structure_data, data) ==
          BelowRadialExponential([0, 0, 1], [0, 0, 0], 5.0, 2.0, 0.5)
end

@testitem "RadialExponential gradients" begin
    using ForwardDiff
    using StaticArrays
    # Away from r = 0 (the cone's apex), the surface is smooth: compare
    # against ForwardDiff.
    points = [
        SVector{3,Float64}(1.0, 2.0, 3.0),
        SVector{3,Float64}(0.0, 0.0, 10.0),
        SVector{3,Float64}(-2.0, 3.0, -1.0),
        SVector{3,Float64}(1e-3, 0.0, 0.0),
    ]
    for c in (
        AboveRadialExponential([0, 0, 1], [0, 0, 0], 5.0, 1.0, 0.5),
        BelowRadialExponential([0, 0, 1], [0, 0, 0], 5.0, 1.0, 0.5),
        AboveRadialExponential([1, 1, 1] ./ sqrt(3.0), [0.2, -0.3, 0.1], 1.0, -2.0, -0.8),
    )
        for x in points
            @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
        end
    end
    # At r = 0 (the cone's genuine kink), check the analytic gradient is
    # finite (not NaN/Inf from a 0/0 division) rather than comparing to
    # ForwardDiff, which is itself singular there.
    c = AboveRadialExponential([0, 0, 1], [0, 0, 0], 5.0, 1.0, 0.5)
    x = SVector{3,Float64}(0.0, 0.0, 0.0)
    g = Packmol.constraint_gradient(c, x)
    @test all(isfinite, g)
end
