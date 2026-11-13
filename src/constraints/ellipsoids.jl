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

"""
    InsideEllipsoid(center, a, b, c, scale; weight)
    InsideEllipsoid(; center, a, b, c, scale, weight)

Keeps the atoms inside an axis-aligned ellipsoid, as the input file line
`inside ellipsoid cx cy cz a b c scale`: semi-axes `a`, `b` and `c`, all multiplied by
`scale`.
"""
InsideEllipsoid(args...; kwargs...) = Ellipsoid{Inside}(args...; kwargs...)
"""
    OutsideEllipsoid(center, a, b, c, scale; weight)
    OutsideEllipsoid(; center, a, b, c, scale, weight)

Keeps the atoms outside an axis-aligned ellipsoid, given as in [`InsideEllipsoid`](@ref).
"""
OutsideEllipsoid(args...; kwargs...) = Ellipsoid{Outside}(args...; kwargs...)

# The violation, about the distance (Å) to the surface: `m (ρ - |scale|)`, for ρ² =
# Σ ((x - center)ᵢ/aᵢ)² and `m` the smallest semi-axis; > 0 outside, < 0 inside, the exact
# distance for a sphere. Not Fortran Packmol's `ρ² - scale²`, which has no units: with
# semi-axes of hundreds of Å, the gradient of a violation of a few hundredths of Å was
# below the stopping tolerance of the optimizer (1e-5) while the fits of the initial
# placement asked for penalties below 1e-10, and stalled (an LNP of 70 nm: 200 loops).
function _ellipsoid_w(c::Ellipsoid, x)
    dx = x - c.center
    ρ = sqrt((dx[1]/c.a)^2 + (dx[2]/c.b)^2 + (dx[3]/c.c)^2)
    return min(c.a, c.b, c.c) * (ρ - abs(c.scale))
end
function _ellipsoid_gradw(c::Ellipsoid, x)
    dx = x - c.center
    q = SVector(dx[1]/c.a^2, dx[2]/c.b^2, dx[3]/c.c^2)
    ρ = sqrt((dx[1]/c.a)^2 + (dx[2]/c.b)^2 + (dx[3]/c.c)^2)
    # (zero at the center, where ρ is not differentiable)
    return ρ > 0 ? min(c.a, c.b, c.c) * q / ρ : zero(q)
end

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
# Input parsing functions: methods of `parse_constraint`, whose keywords must also be
# listed in `constraint_keywords` (constraint_types.jl):
#
function parse_constraint(::Val{Symbol("inside ellipsoid")}, data::AbstractVector{<:AbstractString}, ::Type{T}, filename::String) where {T}
    center, a, b, c, scale = try
        parse.(T, data[1:3]), parse(T, data[4]), parse(T, data[5]), parse(T, data[6]), parse(T, data[7])
    catch
        error("Error parsing 'inside ellipsoid' constraint data for $filename.")
    end
    return Ellipsoid{Inside,T}(; center, a, b, c, scale)
end

function parse_constraint(::Val{Symbol("outside ellipsoid")}, data::AbstractVector{<:AbstractString}, ::Type{T}, filename::String) where {T}
    center, a, b, c, scale = try
        parse.(T, data[1:3]), parse(T, data[4]), parse(T, data[5]), parse(T, data[6]), parse(T, data[7])
    catch
        error("Error parsing 'outside ellipsoid' constraint data for $filename.")
    end
    return Ellipsoid{Outside,T}(; center, a, b, c, scale)
end

# The bounding box of the region (see constraints_base.jl)
function _inside_box(c::Ellipsoid{Inside})
    half = abs(c.scale) * SVector(c.a, c.b, c.c)
    return (c.center - half, c.center + half)
end
