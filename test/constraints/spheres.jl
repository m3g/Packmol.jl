@testitem "Sphere constructors" begin
    # Default weight is weight_default[:sphere] == 1e-2, not Box/Cube/Plane's 5.0.
    @test InsideSphere([0,0,0],1.) == Sphere{Inside,Float64}([0.,0.,0.],1.,1e-2)
    @test InsideSphere(center=[0,0,0],radius=1.) == Sphere{Inside,Float64}([0.,0.,0.],1.,1e-2)
    @test InsideSphere(center=[0,0,0],radius=1.,weight=2.0) == Sphere{Inside,Float64}([0.,0.,0.],1.,2.0)
    @test OutsideSphere([0,0,0],1.) == Sphere{Outside,Float64}([0.,0.,0.],1.,1e-2)
    @test OutsideSphere(center=[0,0,0],radius=1.) == Sphere{Outside,Float64}([0.,0.,0.],1.,1e-2)
    @test OutsideSphere(center=[0,0,0],radius=1.,weight=2.0) == Sphere{Outside,Float64}([0.,0.,0.],1.,2.0)
end

@testitem "Sphere gradients" begin
    using ForwardDiff
    using StaticArrays
    # InsideSphere: point outside the sphere (penalty > 0)
    x = SVector{3,Float64}(1.5, 1.0, 0.)
    c = InsideSphere(center=[0.2, 0., 0.1], radius=0.1)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    # InsideSphere: point inside the sphere (penalty = 0)
    x_in = SVector{3,Float64}(0.21, 0.01, 0.11)
    @test Packmol.constraint_penalty(c, x_in) == 0.0
    @test Packmol.constraint_gradient(c, x_in) == zero(x_in)
    # OutsideSphere: point inside the sphere (penalty > 0)
    c2 = OutsideSphere(center=[0.2, 0., 0.1], radius=1.0)
    x_inside = SVector{3,Float64}(0.3, 0.1, 0.2)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c2, x), x_inside) ≈ Packmol.constraint_gradient(c2, x_inside)
    # OutsideSphere: point outside the sphere (penalty = 0)
    x_out = SVector{3,Float64}(5.0, 5.0, 5.0)
    @test Packmol.constraint_penalty(c2, x_out) == 0.0
    @test Packmol.constraint_gradient(c2, x_out) == zero(x_out)
end
