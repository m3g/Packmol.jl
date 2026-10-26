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

@testitem "_constraint_dispatch" begin
    using StaticArrays
    using Packmol: AnyConstraint, _constraint_dispatch, constraint_penalty, constraint_gradient
    using Packmol: Sphere, Box, Plane, Inside, Outside, Over
    cs = AnyConstraint{Float64}[
        Sphere{Inside,Float64}(center=SVector(0.0, 0.0, 0.0), radius=1.0, weight=5.0),
        Box{Outside,Float64}(center=SVector(0.0, 0.0, 0.0), sides=SVector(1.0, 1.0, 1.0), weight=5.0),
        Plane{Over,Float64}(normal=SVector(0.0, 0.0, 1.0), d=1.0),
    ]
    x = SVector(2.0, 0.1, 0.3)
    for c in cs
        # Called with the concrete type of the constraint
        @test _constraint_dispatch((c, x) -> typeof(c), c, x) === typeof(c)
        @test _constraint_dispatch(constraint_penalty, c, x) == constraint_penalty(c, x)
        @test _constraint_dispatch(constraint_gradient, c, x) == constraint_gradient(c, x)
    end
    f(cs, x) = sum(c -> _constraint_dispatch(constraint_penalty, c, x), cs)
    f(cs, x)
    @test (@allocated f(cs, x)) == 0
end
