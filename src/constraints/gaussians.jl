#
# Gaussian constraints: above/below a Gaussian-shaped surface.
#
# Two flavors, both generalizing `Plane` (they degenerate exactly to it when
# `amplitude == 0`):
#
#   - `Gaussian` (ridge): the surface is `up·x = h(along·x)`, i.e. a Gaussian
#     bump extruded along the direction perpendicular to `along` (within the
#     plane orthogonal to `up`). `up` plays the same role as `Plane`'s
#     `normal` and is used un-normalized, exactly as in `Plane`; `along` is
#     normalized internally (it need not be given as a unit vector).
#   - `RadialGaussian` (dome): the surface is `up·x = h(r)`, where `r` is the
#     radial distance from the axis line through `center` parallel to `up` —
#     a radially symmetric dome rather than an extruded ridge. `up` is
#     normalized internally here, since it doubles as the axis for the
#     radial decomposition (as in `Cylinder`).
#
export Gaussian, AboveGaussian, BelowGaussian
export RadialGaussian, AboveRadialGaussian, BelowRadialGaussian

# Default weights
weight_default[:gaussian] = 5.0
weight_default[:radial_gaussian] = 5.0

#
# Ridge Gaussian
#
@kwdef struct Gaussian{Placement,T} <: Constraint{Placement,3,T}
    up::SVector{3,T}
    along::SVector{3,T}
    d0::T
    amplitude::T
    center::T
    sigma::T
    weight::T = weight_default[:gaussian]
end

function Gaussian{Placement}(up, along, d0, amplitude, center, sigma, weight=weight_default[:gaussian]) where {Placement}
    T = promote_type(eltype(up), eltype(along), typeof(d0), typeof(amplitude), typeof(center), typeof(sigma), typeof(weight))
    return Gaussian{Placement,T}(SVector{3,T}(up), SVector{3,T}(along), T(d0), T(amplitude), T(center), T(sigma), T(weight))
end
Gaussian{Placement}(; up, along, d0, amplitude, center, sigma, weight=weight_default[:gaussian]) where {Placement} =
    Gaussian{Placement}(up, along, d0, amplitude, center, sigma, weight)

AboveGaussian(args...; kargs...) = Gaussian{Over}(args...; kargs...)
BelowGaussian(args...; kargs...) = Gaussian{Below}(args...; kargs...)

# Profile value and derivative at s = along_unit·x
function _gaussian_h(c::Gaussian, s)
    (; d0, amplitude, center, sigma) = c
    return d0 + amplitude * exp(-(s - center)^2 / (2 * sigma^2))
end
_gaussian_dh(c::Gaussian, s, h) = -(s - c.center) / c.sigma^2 * (h - c.d0)

function _gaussian_v_gradv(c::Gaussian, x)
    a = c.along / norm(c.along)
    s = dot(a, x)
    h = _gaussian_h(c, s)
    v = dot(c.up, x) - h
    dh = _gaussian_dh(c, s, h)
    gradv = c.up - dh * a
    return v, gradv
end

_op_gaussian(::Gaussian{Over}, v) = v < zero(v)
_op_gaussian(::Gaussian{Below}, v) = v > zero(v)

function constraint_penalty(c::Gaussian, x)
    v, _ = _gaussian_v_gradv(c, x)
    if _op_gaussian(c, v)
        return c.weight * v^2
    else
        return zero(eltype(x))
    end
end

function constraint_gradient(c::Gaussian, x)
    v, gradv = _gaussian_v_gradv(c, x)
    if _op_gaussian(c, v)
        return 2 * c.weight * v * gradv
    else
        return zero(x)
    end
end

#
# Radial Gaussian (dome)
#
@kwdef struct RadialGaussian{Placement,T} <: Constraint{Placement,3,T}
    up::SVector{3,T}
    center::SVector{3,T}
    d0::T
    amplitude::T
    sigma::T
    weight::T = weight_default[:radial_gaussian]
end

function RadialGaussian{Placement}(up, center, d0, amplitude, sigma, weight=weight_default[:radial_gaussian]) where {Placement}
    T = promote_type(eltype(up), eltype(center), typeof(d0), typeof(amplitude), typeof(sigma), typeof(weight))
    return RadialGaussian{Placement,T}(SVector{3,T}(up), SVector{3,T}(center), T(d0), T(amplitude), T(sigma), T(weight))
end
RadialGaussian{Placement}(; up, center, d0, amplitude, sigma, weight=weight_default[:radial_gaussian]) where {Placement} =
    RadialGaussian{Placement}(up, center, d0, amplitude, sigma, weight)

AboveRadialGaussian(args...; kargs...) = RadialGaussian{Over}(args...; kargs...)
BelowRadialGaussian(args...; kargs...) = RadialGaussian{Below}(args...; kargs...)

# Decompose x - center into the axial coordinate `w` (along the normalized
# `up`) and the perpendicular vector `perp` (so `rsq = sum(abs2, perp)` is
# the squared radial distance to the up-axis line through `center`).
function _radial_gaussian_axial_perp(c::RadialGaussian, x)
    u = c.up / norm(c.up)
    a = x - c.center
    w = dot(a, u)
    perp = a - w * u
    return u, w, perp
end

function _radial_gaussian_v_gradv(c::RadialGaussian, x)
    u, w, perp = _radial_gaussian_axial_perp(c, x)
    rsq = sum(abs2, perp)
    h = c.d0 + c.amplitude * exp(-rsq / (2 * c.sigma^2))
    v = w - h
    dh_drsq = -(h - c.d0) / (2 * c.sigma^2)
    gradv = u - dh_drsq * (2 * perp)
    return v, gradv
end

_op_gaussian(::RadialGaussian{Over}, v) = v < zero(v)
_op_gaussian(::RadialGaussian{Below}, v) = v > zero(v)

function constraint_penalty(c::RadialGaussian, x)
    v, _ = _radial_gaussian_v_gradv(c, x)
    if _op_gaussian(c, v)
        return c.weight * v^2
    else
        return zero(eltype(x))
    end
end

function constraint_gradient(c::RadialGaussian, x)
    v, gradv = _radial_gaussian_v_gradv(c, x)
    if _op_gaussian(c, v)
        return 2 * c.weight * v * gradv
    else
        return zero(x)
    end
end

#
# Input parsing functions: must be appended to the "parse_constraint" dictionary:
#
parse_constraint["above gaussian"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, along, d0, amplitude, center, sigma = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9]), parse(T, data[10])
    catch
        error("Error parsing 'above gaussian' constraint data for $(structure_data[:filename]).")
    end
    return Gaussian{Over,T}(; up, along, d0, amplitude, center, sigma)
end
parse_constraint["below gaussian"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, along, d0, amplitude, center, sigma = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9]), parse(T, data[10])
    catch
        error("Error parsing 'below gaussian' constraint data for $(structure_data[:filename]).")
    end
    return Gaussian{Below,T}(; up, along, d0, amplitude, center, sigma)
end

parse_constraint["above radial_gaussian"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, center, d0, amplitude, sigma = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9])
    catch
        error("Error parsing 'above radial_gaussian' constraint data for $(structure_data[:filename]).")
    end
    return RadialGaussian{Over,T}(; up, center, d0, amplitude, sigma)
end
parse_constraint["below radial_gaussian"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    up, center, d0, amplitude, sigma = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8]), parse(T, data[9])
    catch
        error("Error parsing 'below radial_gaussian' constraint data for $(structure_data[:filename]).")
    end
    return RadialGaussian{Below,T}(; up, center, d0, amplitude, sigma)
end
