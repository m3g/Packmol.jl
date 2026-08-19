#
# `AnyConstraint{T}` is a closed Union of every concrete constraint type,
# used as the element type of `StructureType.constraints` instead of the
# abstract `Constraint`.
#
# This matters for performance: a `Vector{Constraint}` (abstract eltype)
# forces every `constraint_penalty`/`constraint_gradient` call in the hot
# packing loop through dynamic dispatch, which — because the compiler can't
# bound the return type across an open set of possible methods — heap-boxes
# the returned `SVector` on every single call. A `Vector{AnyConstraint{T}}`
# is a small, closed Union of isbits types, which Julia stores inline and
# dispatches via "union splitting" instead: allocation-free, and several
# times faster (measured ~0 vs ~100 bytes/call for a 2-constraint mix).
#
# New constraint types must be added to this Union as well as to
# `parse_constraint`/exports in their own file.
#
const AnyConstraint{T} = Union{
    Box{Inside,T}, Box{Outside,T},
    Cube{Inside,T}, Cube{Outside,T},
    Sphere{Inside,T}, Sphere{Outside,T},
    Plane{Over,T}, Plane{Below,T},
    Cylinder{Inside,T}, Cylinder{Outside,T},
    Ellipsoid{Inside,T}, Ellipsoid{Outside,T},
}

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
