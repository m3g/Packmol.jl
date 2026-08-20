@testitem "Ellipsoid constructors" begin
    @test InsideEllipsoid([0,0,0],2.,1.,1.,1.) == Ellipsoid{Inside,Float64}([0.,0.,0.],2.,1.,1.,1.,5.0)
    @test InsideEllipsoid(center=[0,0,0],a=2.,b=1.,c=1.,scale=1.) == Ellipsoid{Inside,Float64}([0.,0.,0.],2.,1.,1.,1.,5.0)
    @test InsideEllipsoid(center=[0,0,0],a=2.,b=1.,c=1.,scale=1.,weight=2.0) == Ellipsoid{Inside,Float64}([0.,0.,0.],2.,1.,1.,1.,2.0)
    @test OutsideEllipsoid([0,0,0],2.,1.,1.,1.) == Ellipsoid{Outside,Float64}([0.,0.,0.],2.,1.,1.,1.,5.0)
    @test OutsideEllipsoid(center=[0,0,0],a=2.,b=1.,c=1.,scale=1.) == Ellipsoid{Outside,Float64}([0.,0.,0.],2.,1.,1.,1.,5.0)
end

@testitem "Ellipsoid gradients" begin
    using ForwardDiff
    using StaticArrays
    c = InsideEllipsoid(center=[0.2, -0.3, 0.1], a=2.0, b=1.0, c=1.5, scale=1.0)
    for x in (
        SVector(0.3, 0.1, 0.2),   # inside
        SVector(3.0, 0.1, 0.2),   # outside along a
        SVector(0.3, 2.0, 0.2),   # outside along b
        SVector(0.3, 0.1, 3.0),   # outside along c
    )
        @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    end

    c2 = OutsideEllipsoid(center=[0.2, -0.3, 0.1], a=2.0, b=1.0, c=1.5, scale=1.0)
    for x in (
        SVector(0.3, 0.1, 0.2),   # inside (forbidden)
        SVector(3.0, 0.1, 0.2),   # outside (already fine)
    )
        @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c2, x), x) ≈ Packmol.constraint_gradient(c2, x)
    end

    # scale != 1: effective semi-axes are (a*scale, b*scale, c*scale)
    c3 = InsideEllipsoid(center=[0.0, 0.0, 0.0], a=1.0, b=1.0, c=1.0, scale=2.0)
    @test Packmol.constraint_penalty(c3, SVector(1.5, 0.0, 0.0)) == 0.0
    @test Packmol.constraint_penalty(c3, SVector(2.5, 0.0, 0.0)) > 0.0
    for x in (SVector(1.5, 0.0, 0.0), SVector(2.5, 0.0, 0.0))
        @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c3, x), x) ≈ Packmol.constraint_gradient(c3, x)
    end
end
