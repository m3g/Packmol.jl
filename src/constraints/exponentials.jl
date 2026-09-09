#
# Exponential constraints: above/below an exponential surface.
#
# Two flavors, both generalizing `Plane` (they degenerate exactly to it when
# `amplitude == 0`), mirroring `Gaussian`/`RadialGaussian` in gaussians.jl:
#
#   - `Exponential` (ridge): the surface is `up·x = h(along·x)`, an
#     exponential ramp extruded along the direction perpendicular to `along`
#     (within the plane orthogonal to `up`). `up` is used un-normalized, as
#     in `Plane`; `along` is normalized internally. Smooth everywhere.
#   - `RadialExponential` (cone/spike): the surface is `up·x = h(r)`, where
#     `r` is the radial distance from the axis line through `center`
#     parallel to `up`. `up` is normalized internally, since it doubles as
#     the axis for the radial decomposition (as in `Cylinder`). `h` is
#     evaluated on `r` itself (not `r^2`), which is what gives this variant
#     its characteristic cone-like point (a genuine kink, not just a
#     removable singularity) at `r == 0`; the gradient there is taken to be
#     purely axial, as for `RadialWave`.
#
export Exponential, AboveExponential, BelowExponential
export RadialExponential, AboveRadialExponential, BelowRadialExponential

# Default weights
weight_default[:exponential] = 5.0
weight_default[:radial_exponential] = 5.0

#
# Ridge exponential
#
@kwdef struct Exponential{Placement,T} <: Constraint{Placement,3,T}
    up::SVector{3,T}
    along::SVector{3,T}
    d0::T
    amplitude::T
    center::T
    rate::T
    weight::T = weight_default[:exponential]
end

function Exponential{Placement}(up, along, d0, amplitude, center, rate, weight=weight_default[:exponential]) where {Placement}
    T = promote_type(eltype(up), eltype(along), typeof(d0), typeof(amplitude), typeof(center), typeof(rate), typeof(weight))
    return Exponential{Placement,T}(SVector{3,T}(up), SVector{3,T}(along), T(d0), T(amplitude), T(center), T(rate), T(weight))
end
Exponential{Placement}(; up, along, d0, amplitude, center, rate, weight=weight_default[:exponential]) where {Placement} =
    Exponential{Placement}(up, along, d0, amplitude, center, rate, weight)

AboveExponential(args...; kargs...) = Exponential{Over}(args...; kargs...)
BelowExponential(args...; kargs...) = Exponential{Below}(args...; kargs...)

# Profile value and derivative at s = along_unit·x
function _exponential_h(c::Exponential, s)
    (; d0, amplitude, center, rate) = c
    return d0 + amplitude * exp(rate * (s - center))
end
_exponential_dh(c::Exponential, h) = c.rate * (h - c.d0)

function _exponential_v_gradv(c::Exponential, x)
    a = c.along / norm(c.along)
    s = dot(a, x)
    h = _exponential_h(c, s)
    v = dot(c.up, x) - h
    dh = _exponential_dh(c, h)
    gradv = c.up - dh * a
    return v, gradv
end

_op_exponential(::Exponential{Over}, v) = v < zero(v)
_op_exponential(::Exponential{Below}, v) = v > zero(v)

function constraint_penalty(c::Exponential, x)
    v, _ = _exponential_v_gradv(c, x)
    if _op_exponential(c, v)
        return c.weight * v^2
    else
        return zero(eltype(x))
    end
end

function constraint_gradient(c::Exponential, x)
    v, gradv = _exponential_v_gradv(c, x)
    if _op_exponential(c, v)
        return 2 * c.weight * v * gradv
    else
        return zero(x)
    end
end

#
# Radial exponential (cone/spike)
#
@kwdef struct RadialExponential{Placement,T} <: Constraint{Placement,3,T}
    up::SVector{3,T}
    center::SVector{3,T}
    d0::T
    amplitude::T
    rate::T
    weight::T = weight_default[:radial_exponential]
end

function RadialExponential{Placement}(up, center, d0, amplitude, rate, weight=weight_default[:radial_exponential]) where {Placement}
    T = promote_type(eltype(up), eltype(center), typeof(d0), typeof(amplitude), typeof(rate), typeof(weight))
    return RadialExponential{Placement,T}(SVector{3,T}(up), SVector{3,T}(center), T(d0), T(amplitude), T(rate), T(weight))
end
RadialExponential{Placement}(; up, center, d0, amplitude, rate, weight=weight_default[:radial_exponential]) where {Placement} =
    RadialExponential{Placement}(up, center, d0, amplitude, rate, weight)

AboveRadialExponential(args...; kargs...) = RadialExponential{Over}(args...; kargs...)
BelowRadialExponential(args...; kargs...) = RadialExponential{Below}(args...; kargs...)

# Decompose x - center into the axial coordinate `w` (along the normalized
# `up`) and the perpendicular vector `perp`, whose norm `r` is the radial
# distance to the up-axis line through `center`.
function _radial_exponential_axial_perp(c::RadialExponential, x)
    u = c.up / norm(c.up)
    a = x - c.center
    w = dot(a, u)
    perp = a - w * u
    return u, w, perp
end

function _radial_exponential_v_gradv(c::RadialExponential, x)
    u, w, perp = _radial_exponential_axial_perp(c, x)
    r = norm(perp)
    h = c.d0 + c.amplitude * exp(c.rate * r)
    v = w - h
    dh_dr = c.rate * (h - c.d0)
    # Radial direction is undefined exactly at r == 0 (the cone's apex);
    # treat its contribution to the gradient as zero there.
    rdir = iszero(r) ? zero(perp) : perp / r
    gradv = u - dh_dr * rdir
    return v, gradv
end

_op_exponential(::RadialExponential{Over}, v) = v < zero(v)
_op_exponential(::RadialExponential{Below}, v) = v > zero(v)

function constraint_penalty(c::RadialExponential, x)
    v, _ = _radial_exponential_v_gradv(c, x)
    if _op_exponential(c, v)
        return c.weight * v^2
    else
        return zero(eltype(x))
    end
end

function constraint_gradient(c::RadialExponential, x)
    v, gradv = _radial_exponential_v_gradv(c, x)
    if _op_exponential(c, v)
        return 2 * c.weight * v * gradv
    else
        return zero(x)
    end
end

#
# Input parsing functions: must be appended to the "parse_constraint" dictionary:
#
parse_constraint["above exponential"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, along, d0, amplitude, center, rate = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9]), parse(T, data[10])
    catch
        error("Error parsing 'above exponential' constraint data for $(structure_data[:filename]).")
    end
    return Exponential{Over,T}(; up, along, d0, amplitude, center, rate)
end
parse_constraint["below exponential"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, along, d0, amplitude, center, rate = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9]), parse(T, data[10])
    catch
        error("Error parsing 'below exponential' constraint data for $(structure_data[:filename]).")
    end
    return Exponential{Below,T}(; up, along, d0, amplitude, center, rate)
end

parse_constraint["above radial_exponential"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, center, d0, amplitude, rate = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9])
    catch
        error("Error parsing 'above radial_exponential' constraint data for $(structure_data[:filename]).")
    end
    return RadialExponential{Over,T}(; up, center, d0, amplitude, rate)
end
parse_constraint["below radial_exponential"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, center, d0, amplitude, rate = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9])
    catch
        error("Error parsing 'below radial_exponential' constraint data for $(structure_data[:filename]).")
    end
    return RadialExponential{Below,T}(; up, center, d0, amplitude, rate)
end
