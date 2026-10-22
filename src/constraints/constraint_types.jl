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
# New constraint types must be added to this Union, their input keywords to
# `constraint_keywords` below, and `parse_constraint` methods/exports to their
# own file.
#
const AnyConstraint{T} = Union{
    Box{Inside,T}, Box{Outside,T},
    Cube{Inside,T}, Cube{Outside,T},
    Sphere{Inside,T}, Sphere{Outside,T},
    Plane{Over,T}, Plane{Below,T},
    Cylinder{Inside,T}, Cylinder{Outside,T},
    Ellipsoid{Inside,T}, Ellipsoid{Outside,T},
    Gaussian{Over,T}, Gaussian{Below,T},
    RadialGaussian{Over,T}, RadialGaussian{Below,T},
    Wave{Over,T}, Wave{Below,T},
    RadialWave{Over,T}, RadialWave{Below,T},
    Exponential{Over,T}, Exponential{Below,T},
    RadialExponential{Over,T}, RadialExponential{Below,T},
}

#
# `_constraint_dispatch(f, c, args...)` calls `f(c, args...)` for a constraint `c` of
# the `AnyConstraint` union through an explicit chain of `isa` branches, one per type of
# the union, such that `f` is always called with the concrete type of `c`. Julia's
# automatic union splitting covers only small unions (up to 4 types), and a dynamic
# dispatch does not compile with `juliac --trim` (and allocates). Every call with a
# constraint taken from a `Vector{AnyConstraint{T}}` must go through this function.
#
let branches = nothing
    for t in reverse(Base.uniontypes(AnyConstraint{Float64}))
        concrete = Expr(:curly, t.name.wrapper, t.parameters[1], :T)
        call = :(f(c::$concrete, args...))
        branches = isnothing(branches) ? call : Expr(:if, :(c isa $concrete), call, branches)
    end
    @eval @inline function _constraint_dispatch(f::F, c::AnyConstraint{T}, args...) where {F,T}
        $branches
    end
end

#
# Input file keywords of the constraints, each parsed by its own
# `parse_constraint(::Val{Symbol(keyword)}, ...)` method.
#
const constraint_keywords = (
    "inside box", "outside box", "inside cube", "outside cube",
    "inside sphere", "outside sphere",
    "above plane", "over plane", "below plane",
    "inside cylinder", "outside cylinder",
    "inside ellipsoid", "outside ellipsoid",
    "above gaussian", "below gaussian", "above radial_gaussian", "below radial_gaussian",
    "above sin", "below sin", "above cos", "below cos",
    "above radial_sin", "below radial_sin", "above radial_cos", "below radial_cos",
    "above exponential", "below exponential", "above radial_exponential", "below radial_exponential",
)
const _constraint_vals = map(keyword -> Val(Symbol(keyword)), constraint_keywords)

# The first word of the constraint keywords (e.g. "inside")
const constraint_placements = unique(String[first(split(keyword)) for keyword in constraint_keywords])

#
# Parse the data of the constraint given by `keyword` (e.g. "inside box") and
# push it to `constraints`. A recursion over the keywords, rather than a lookup,
# so that each constraint is built and pushed with its own concrete type.
#
function push_parsed_constraint!(
    constraints::Vector{AnyConstraint{T}}, keyword::AbstractString, data, filename::String,
    vals::Tuple=_constraint_vals,
) where {T}
    isempty(vals) && throw(ArgumentError("Constraint $keyword not recognized, in structure $filename."))
    if keyword == String(_unval(first(vals)))
        push!(constraints, parse_constraint(first(vals), data, T, filename))
    else
        push_parsed_constraint!(constraints, keyword, data, filename, Base.tail(vals))
    end
    return constraints
end
_unval(::Val{S}) where {S} = S
