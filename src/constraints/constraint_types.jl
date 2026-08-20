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
