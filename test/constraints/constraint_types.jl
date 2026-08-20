@testitem "AnyConstraint covers every constraint type" begin
    using StaticArrays
    T = Float64
    constraints = Packmol.AnyConstraint{T}[
        InsideBox([0.,0.,0.],[1.,1.,1.]),
        OutsideBox([0.,0.,0.],[1.,1.,1.]),
        InsideCube([0.,0.,0.],1.),
        OutsideCube([0.,0.,0.],1.),
        InsideSphere([0.,0.,0.],1.),
        OutsideSphere([0.,0.,0.],1.),
        AbovePlane([0.,0.,1.],1.),
        BelowPlane([0.,0.,1.],1.),
        InsideCylinder([0.,0.,0.],[0.,0.,1.],1.,1.),
        OutsideCylinder([0.,0.,0.],[0.,0.,1.],1.,1.),
        InsideEllipsoid([0.,0.,0.],1.,1.,1.,1.),
        OutsideEllipsoid([0.,0.,0.],1.,1.,1.,1.),
    ]
    @test eltype(constraints) == Packmol.AnyConstraint{T}
    x = SVector(0.5,0.5,0.5)
    # No error / correct dispatch for every stored type.
    @test all(c -> Packmol.constraint_penalty(c, x) isa T, constraints)
    @test all(c -> Packmol.constraint_gradient(c, x) isa SVector{3,T}, constraints)
end
