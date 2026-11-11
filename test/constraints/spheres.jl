@testitem "Sphere constructors" begin
    # Default weight 500.0 (the violation is a distance in Å; see spheres.jl)
    @test InsideSphere([0,0,0],1.) == Sphere{Inside,Float64}([0.,0.,0.],1.,500.0)
    @test InsideSphere(center=[0,0,0],radius=1.) == Sphere{Inside,Float64}([0.,0.,0.],1.,500.0)
    @test InsideSphere(center=[0,0,0],radius=1.,weight=2.0) == Sphere{Inside,Float64}([0.,0.,0.],1.,2.0)
    @test OutsideSphere([0,0,0],1.) == Sphere{Outside,Float64}([0.,0.,0.],1.,500.0)
    @test OutsideSphere(center=[0,0,0],radius=1.) == Sphere{Outside,Float64}([0.,0.,0.],1.,500.0)
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

@testitem "Sphere penalty scale" begin
    # The violation is the distance to the surface (Å): penalty weight·δ², gradient
    # 2·weight·δ along the radius, whatever the radius (it was d² - r², whose square
    # grows as d⁴: 10¹⁶ for atoms hundreds of Å out of a sphere of 890 Å)
    using StaticArrays
    using LinearAlgebra: norm
    for R in (10.0, 890.0)
        c = InsideSphere(center = [1.0, 2.0, 3.0], radius = R, weight = 5.0)
        x = SVector(1.0, 2.0, 3.0 + R + 0.5)
        @test Packmol.constraint_penalty(c, x) ≈ 5.0 * 0.25
        @test Packmol.constraint_gradient(c, x) ≈ SVector(0.0, 0.0, 5.0)
        @test Packmol.constraint_penalty(c, SVector(1.0, 2.0, 3.0 + 3R)) ≈ 5.0 * (2R)^2
        o = OutsideSphere(center = [1.0, 2.0, 3.0], radius = R, weight = 5.0)
        @test Packmol.constraint_penalty(o, SVector(1.0, 2.0, 3.0 + R - 0.5)) ≈ 5.0 * 0.25
        @test Packmol.constraint_gradient(o, SVector(1.0, 2.0, 3.0 + R - 0.5)) ≈ SVector(0.0, 0.0, -5.0)
        # at the center: penalty weight·R², zero gradient (not differentiable)
        @test Packmol.constraint_penalty(o, SVector(1.0, 2.0, 3.0)) ≈ 5.0 * R^2
        @test Packmol.constraint_gradient(o, SVector(1.0, 2.0, 3.0)) == zero(SVector{3,Float64})
    end
end
