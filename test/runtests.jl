using Packmol
using Test
using ForwardDiff
using StaticArrays

@testset "Constraints" begin

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
