#
# Spherical constraints
#
export Sphere, InsideSphere, OutsideSphere

# Default weights
weight_default[:sphere] = 500.0

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

"""
    InsideSphere(center, radius; weight)
    InsideSphere(; center, radius, weight)

Keeps the atoms inside a sphere, as the input file line `inside sphere cx cy cz radius`.
"""
InsideSphere(args...; kwargs...) = Sphere{Inside}(args...; kwargs...)
"""
    OutsideSphere(center, radius; weight)
    OutsideSphere(; center, radius, weight)

Keeps the atoms outside a sphere, given as in [`InsideSphere`](@ref).
"""
OutsideSphere(args...; kwargs...) = Sphere{Outside}(args...; kwargs...)

# The violation is the distance (Å) to the surface, `d - radius` (> 0 outside, < 0
# inside), and the penalty `weight * violation^2`. Not Fortran Packmol's `d² - radius²`,
# in Å²: its square grows as d⁴, and atoms of a long molecule placed hundreds of Å out of
# a large sphere (a cell of 890 Å) gave an objective of 10¹⁶ and a first step of the
# optimizer that threw them to 10¹⁰ Å. Near the surface, `d² - r²` ≈ 2r (d - r): the old
# weight (0.01) gave a stiffness of 0.04 r², which the default weight (500) matches at
# r ≈ 110 Å. A stiff wall matters for packings of many molecules against a sphere: with
# a weight of 5, the exosome needed 149 loops (59 with 500), and with 5000 Zika needed 21
# (16 with 500).
_sphere_w(c::Sphere, x) = norm(x - c.center) - c.radius
# Gradient of the distance: the radial direction (zero at the center, where it is not
# differentiable)
function _sphere_gradw(c::Sphere, x)
    dx = x - c.center
    d = norm(dx)
    return d > 0 ? dx / d : zero(dx)
end

function constraint_penalty(c::Sphere{Inside}, x)
    w = _sphere_w(c, x)
    return w > 0 ? c.weight * w^2 : zero(eltype(x))
end

function constraint_penalty(c::Sphere{Outside}, x)
    w = _sphere_w(c, x)
    return w < 0 ? c.weight * w^2 : zero(eltype(x))
end

function constraint_gradient(c::Sphere{Inside}, x)
    w = _sphere_w(c, x)
    return w > 0 ? 2 * c.weight * w * _sphere_gradw(c, x) : zero(x)
end

function constraint_gradient(c::Sphere{Outside}, x)
    w = _sphere_w(c, x)
    return w < 0 ? 2 * c.weight * w * _sphere_gradw(c, x) : zero(x)
end

#
# Input parsing functions: methods of `parse_constraint`, whose keywords must also be
# listed in `constraint_keywords` (constraint_types.jl):
#
function parse_constraint(::Val{Symbol("inside sphere")}, data::AbstractVector{<:AbstractString}, ::Type{T}, filename::String) where {T}
    center, radius = try
        parse.(T, data[1:3]), parse.(T, data[4])
    catch
        error("Error parsing 'inside sphere' constraint data for $filename.")
    end
    return Sphere{Inside,T}(;center, radius)
end

function parse_constraint(::Val{Symbol("outside sphere")}, data::AbstractVector{<:AbstractString}, ::Type{T}, filename::String) where {T}
    center, radius = try
        parse.(T, data[1:3]), parse.(T, data[4])
    catch
        error("Error parsing 'outside sphere' constraint data for $filename.")
    end
    return Sphere{Outside,T}(;center, radius)
end

# The bounding box of the region (see constraints_base.jl)
_inside_box(c::Sphere{Inside}) = (c.center .- c.radius, c.center .+ c.radius)
