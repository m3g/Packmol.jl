using Packmol
using Test
using ForwardDiff
using StaticArrays

@testset "Constraint constructors" begin

    @test InsideCube([0,0,0],1.) == Cube{Inside,3,Float64}([0.,0.,0.],1.,5.0)
    @test InsideCube(center=[0,0,0],side=1.) == Cube{Inside,3,Float64}([0.,0.,0.],1.,5.0)
    @test InsideCube(center=[0,0,0],side=1.,weight=2.0) == Cube{Inside,3,Float64}([0.,0.,0.],1.,2.0)
    @test OutsideCube([0,0,0],1.) == Cube{Outside,3,Float64}([0.,0.,0.],1.,5.0)
    @test OutsideCube(center=[0,0,0],side=1.) == Cube{Outside,3,Float64}([0.,0.,0.],1.,5.0)
    @test OutsideCube(center=[0,0,0],side=1.,weight=2.0) == Cube{Outside,3,Float64}([0.,0.,0.],1.,2.0)

    @test InsideBox([0,0,0],[1.,1.,1.]) == Box{Inside,3,Float64}([0.,0.,0.],[1.,1.,1.],5.0)
    @test InsideBox(center=[0,0,0],sides=[1.,1.,1.]) == Box{Inside,3,Float64}([0.,0.,0.],[1.,1.,1.],5.0)
    @test InsideBox(center=[0,0,0],sides=[1.,1.,1.],weight=2.0) == Box{Inside,3,Float64}([0.,0.,0.],[1.,1.,1.],2.0)
    @test OutsideBox([0,0,0],[1.,1.,1.]) == Box{Outside,3,Float64}([0.,0.,0.],[1.,1.,1.],5.0)
    @test OutsideBox(center=[0,0,0],sides=[1.,1.,1.]) == Box{Outside,3,Float64}([0.,0.,0.],[1.,1.,1.],5.0)
    @test OutsideBox(center=[0,0,0],sides=[1.,1.,1.],weight=2.0) == Box{Outside,3,Float64}([0.,0.,0.],[1.,1.,1.],2.0)

    @test InsideSphere([0,0,0],1.) == Sphere{Inside,3,Float64}([0.,0.,0.],1.,5.0)
    @test InsideSphere(center=[0,0,0],radius=1.) == Sphere{Inside,3,Float64}([0.,0.,0.],1.,5.0)
    @test InsideSphere(center=[0,0,0],radius=1.,weight=2.0) == Sphere{Inside,3,Float64}([0.,0.,0.],1.,2.0)
    @test OutsideSphere([0,0,0],1.) == Sphere{Outside,3,Float64}([0.,0.,0.],1.,5.0)
    @test OutsideSphere(center=[0,0,0],radius=1.) == Sphere{Outside,3,Float64}([0.,0.,0.],1.,5.0)
    @test OutsideSphere(center=[0,0,0],radius=1.,weight=2.0) == Sphere{Outside,3,Float64}([0.,0.,0.],1.,2.0)

end

@testset "Constraint gradients" begin

    x = SVector{3,Float64}(1.5,1.0,0.)
    c = InsideCube([0.2, 0., 0.1], 0.5)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
    c = OutsideCube([0.2, 0., 0.1], 2.)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
    c = InsideBox([0.2, 0., 0.1], [0.5, 0.7, 1.0])
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
    c = OutsideBox([0.2, 0., 0.1], [0.5, 0.7, 1.0])
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
    c = InsideSphere([0.2, 0., 0.1], 0.1)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
    c = OutsideSphere([0.2, 0., 0.1], 1.0)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)

end
