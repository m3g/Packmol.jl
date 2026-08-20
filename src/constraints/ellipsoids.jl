#
# Ellipsoid constraints: an axis-aligned ellipsoid defined by its `center`,
# semi-axes (`a`, `b`, `c`), and a dimensionless `scale` factor applied
# uniformly to all three semi-axes (so the effective semi-axes are
# `a*scale`, `b*scale`, `c*scale`).
#
# Matches the original Fortran Packmol's `inside/outside ellipsoid` keyword:
#   inside ellipsoid  cx cy cz  a b c  scale
#   outside ellipsoid cx cy cz  a b c  scale
#
export Ellipsoid, InsideEllipsoid, OutsideEllipsoid

# Default weight
weight_default[:ellipsoid] = 5.0

@kwdef struct Ellipsoid{Placement,T} <: Constraint{Placement,3,T}
    center::SVector{3,T}
    a::T
    b::T
    c::T
    scale::T
    weight::T = weight_default[:ellipsoid]
end

# Outer constructors that infer T from the given arguments, so that
# `Ellipsoid{Placement}(...)` (a partially-applied UnionAll, as used by
# InsideEllipsoid/OutsideEllipsoid below) can be called without the type
# parameter T needing to be given explicitly.
function Ellipsoid{Placement}(center, a, b, c, scale, weight=weight_default[:ellipsoid]) where {Placement}
    T = promote_type(eltype(center), typeof(a), typeof(b), typeof(c), typeof(scale), typeof(weight))
    return Ellipsoid{Placement,T}(SVector{3,T}(center), T(a), T(b), T(c), T(scale), T(weight))
end
Ellipsoid{Placement}(; center, a, b, c, scale, weight=weight_default[:ellipsoid]) where {Placement} =
    Ellipsoid{Placement}(center, a, b, c, scale, weight)

InsideEllipsoid(args...; kargs...) = Ellipsoid{Inside}(args...; kargs...)
OutsideEllipsoid(args...; kargs...) = Ellipsoid{Outside}(args...; kargs...)

# Normalized "radius" of x on the ellipsoid's scale: > 0 outside, < 0 inside,
# 0 exactly on the surface `(x-center)/(a,b,c) scaled by `scale`.
function _ellipsoid_w(c::Ellipsoid, x)
    dx = x - c.center
    return (dx[1]/c.a)^2 + (dx[2]/c.b)^2 + (dx[3]/c.c)^2 - c.scale^2
end
_ellipsoid_gradw(c::Ellipsoid, x) = (dx = x - c.center; 2 * SVector(dx[1]/c.a^2, dx[2]/c.b^2, dx[3]/c.c^2))

function constraint_penalty(c::Ellipsoid{Inside}, x)
    w = _ellipsoid_w(c, x)
    if w > 0
        return c.weight * w^2
    else
        return zero(eltype(x))
    end
end

function constraint_penalty(c::Ellipsoid{Outside}, x)
    w = _ellipsoid_w(c, x)
    if w < 0
        return c.weight * w^2
    else
        return zero(eltype(x))
    end
end

function constraint_gradient(c::Ellipsoid{Inside}, x)
    w = _ellipsoid_w(c, x)
    if w > 0
        return 2 * c.weight * w * _ellipsoid_gradw(c, x)
    else
        return zero(x)
    end
end

function constraint_gradient(c::Ellipsoid{Outside}, x)
    w = _ellipsoid_w(c, x)
    if w < 0
        return 2 * c.weight * w * _ellipsoid_gradw(c, x)
    else
        return zero(x)
    end
end

#
# Input parsing functions: must be appended to the "parse_constraint" dictionary:
#
parse_constraint["inside ellipsoid"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    center, a, b, c, scale = try
        parse.(T, data[1:3]), parse(T, data[4]), parse(T, data[5]), parse(T, data[6]), parse(T, data[7])
    catch
        error("Error parsing 'inside ellipsoid' constraint data for $(structure_data[:filename]).")
    end
    return Ellipsoid{Inside,T}(; center, a, b, c, scale)
end

parse_constraint["outside ellipsoid"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    center, a, b, c, scale = try
        parse.(T, data[1:3]), parse(T, data[4]), parse(T, data[5]), parse(T, data[6]), parse(T, data[7])
    catch
        error("Error parsing 'outside ellipsoid' constraint data for $(structure_data[:filename]).")
    end
    return Ellipsoid{Outside,T}(; center, a, b, c, scale)
end
