#
# Spherical constraints
#
export Sphere, InsideSphere, OutsideSphere

# Default weights
weight_default[:sphere] = 1e-2

#
# Spheres
#
@kwdef struct Sphere{Placement,T} <: Constraint{Placement,3,T}
    center::SVector{3,T}
    radius::T
    weight::T = weight_default[:sphere]
end

# Outer constructors that infer T from the given arguments, so that
# `Sphere{Placement}(...)` (a partially-applied UnionAll, as used by
# InsideSphere/OutsideSphere below) can be called without the type
# parameter T needing to be given explicitly.
function Sphere{Placement}(center, radius, weight=weight_default[:sphere]) where {Placement}
    T = promote_type(eltype(center), typeof(radius), typeof(weight))
    return Sphere{Placement,T}(SVector{3,T}(center), T(radius), T(weight))
end
Sphere{Placement}(; center, radius, weight=weight_default[:sphere]) where {Placement} =
    Sphere{Placement}(center, radius, weight)

InsideSphere(args...; kargs...) = Sphere{Inside}(args...; kargs...)
OutsideSphere(args...; kargs...) = Sphere{Outside}(args...; kargs...)

function constraint_penalty(c::Sphere{Inside}, x)
    (; center, radius, weight) = c
    w = sum(abs2, x - center) - radius^2
    if w > 0
        return weight * w^2
    else
        return zero(eltype(x))
    end
end

function constraint_penalty(c::Sphere{Outside}, x)
    (; center, radius, weight) = c
    w = sum(abs2, x - center) - radius^2
    if w < 0
        return weight * w^2
    else
        return zero(eltype(x))
    end
end

function constraint_gradient(c::Sphere{Inside}, x)
    (; center, radius, weight) = c
    dx = x - center
    w = sum(abs2, dx) - radius^2
    if w > 0
        return 4 * weight * w * dx
    else
        return zero(x)
    end
end

function constraint_gradient(c::Sphere{Outside}, x)
    (; center, radius, weight) = c
    dx = x - center
    w = sum(abs2, dx) - radius^2
    if w < 0
        return 4 * weight * w * dx
    else
        return zero(x)
    end
end

#
# Input parsing functions: must be appended to the "parse_constraint" dictionary:
#
parse_constraint["inside sphere"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    center, radius = try
        parse.(T, data[1:3]), parse.(T, data[4])
    catch
        error("Error parsing 'inside sphere' constraint data for $(structure_data[:filename]).")
    end
    return Sphere{Inside,T}(;center, radius)
end

parse_constraint["outside sphere"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    center, radius = try
        parse.(T, data[1:3]), parse.(T, data[4])
    catch
        error("Error parsing 'outside sphere' constraint data for $(structure_data[:filename]).")
    end
    return Sphere{Outside,T}(;center, radius)
end
