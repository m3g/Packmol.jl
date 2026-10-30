#
# Cylinder constraints: a finite, capped cylinder defined by one end of its
# axis (`center`), an axis direction (`axis`, need not be normalized), a
# `radius`, and a `length` along the axis from `center`.
#
# Matches the original Fortran Packmol's `inside/outside cylinder` keyword:
#   inside cylinder  cx cy cz  vx vy vz  radius length
#   outside cylinder cx cy cz  vx vy vz  radius length
#
export Cylinder, InsideCylinder, OutsideCylinder

# Default weight
weight_default[:cylinder] = 5.0

@kwdef struct Cylinder{Placement,T} <: Constraint{Placement,3,T}
    center::SVector{3,T}
    axis::SVector{3,T}
    radius::T
    length::T
    weight::T = weight_default[:cylinder]
end

# Outer constructors that infer T from the given arguments, so that
# `Cylinder{Placement}(...)` (a partially-applied UnionAll, as used by
# InsideCylinder/OutsideCylinder below) can be called without the type
# parameter T needing to be given explicitly.
function Cylinder{Placement}(center, axis, radius, length, weight=weight_default[:cylinder]) where {Placement}
    T = promote_type(eltype(center), eltype(axis), typeof(radius), typeof(length), typeof(weight))
    return Cylinder{Placement,T}(SVector{3,T}(center), SVector{3,T}(axis), T(radius), T(length), T(weight))
end
Cylinder{Placement}(; center, axis, radius, length, weight=weight_default[:cylinder]) where {Placement} =
    Cylinder{Placement}(center, axis, radius, length, weight)

InsideCylinder(args...; kargs...) = Cylinder{Inside}(args...; kargs...)
OutsideCylinder(args...; kargs...) = Cylinder{Outside}(args...; kargs...)

# Decompose x - center into its axial coordinate `w` (signed distance from
# `center` along the normalized axis) and its perpendicular vector `r` (so
# that `sum(abs2, r)` is the squared radial distance from the cylinder axis).
# `r` is, by construction, orthogonal to the axis, which simplifies the
# gradient of `sum(abs2, r)` to `2r` (see below).
function _cylinder_axial_radial(c::Cylinder, x)
    v = c.axis / norm(c.axis)
    a = x - c.center
    w = dot(a, v)
    r = a - w * v
    return v, w, r
end

# Inside: penalty grows if the point is behind the near cap (w < 0), beyond
# the far cap (w > length), or outside the radius — matching the original
# Fortran Packmol's `inside cylinder` penalty (sum of the three violations).
function constraint_penalty(c::Cylinder{Inside}, x)
    _, w, r = _cylinder_axial_radial(c, x)
    t1 = max(-w, zero(w))
    t2 = max(w - c.length, zero(w))
    t3 = max(sum(abs2, r) - c.radius^2, zero(w))
    return c.weight * (t1^2 + t2^2 + t3^2)
end

function constraint_gradient(c::Cylinder{Inside}, x)
    v, w, r = _cylinder_axial_radial(c, x)
    t1 = max(-w, zero(w))
    t2 = max(w - c.length, zero(w))
    t3 = max(sum(abs2, r) - c.radius^2, zero(w))
    return c.weight * (2 * v * (t2 - t1) + 4 * t3 * r)
end

# Outside: zero unless the point is inside the forbidden cylinder (past the near cap,
# before the far cap, and within the radius); then the squared distance to the nearest
# way out (through either cap, or through the side), with the gradient along it. (Fortran
# Packmol uses the product of the three violations, whose scale grows with the length of
# the cylinder: for long cylinders, such as those crossing a periodic box, it is huge in
# the middle of the cylinder, and the packing converges slowly.)
function _outside_cylinder_depth(c::Cylinder, x)
    v, w, r = _cylinder_axial_radial(c, x)
    rr = norm(r)
    depths = (w, c.length - w, c.radius - rr)
    k = argmin(depths)
    depth = depths[k]
    return depth > 0 ? (depth, k) : (zero(depth), 0), v, r, rr
end

function constraint_penalty(c::Cylinder{Outside}, x)
    (depth, _), _ = _outside_cylinder_depth(c, x)
    return c.weight * depth^2
end

function constraint_gradient(c::Cylinder{Outside}, x)
    (depth, k), v, r, rr = _outside_cylinder_depth(c, x)
    k == 0 && return zero(x)
    # gradient of the depth: v through the near cap, -v through the far cap, and
    # minus the radial direction through the side (zero exactly on the axis)
    ddepth = k == 1 ? v : k == 2 ? -v : (rr > 0 ? -r / rr : zero(r))
    return 2 * c.weight * depth * ddepth
end

#
# Input parsing functions: methods of `parse_constraint`, whose keywords must also be
# listed in `constraint_keywords` (constraint_types.jl):
#
function parse_constraint(::Val{Symbol("inside cylinder")}, data::AbstractVector{<:AbstractString}, ::Type{T}, filename::String) where {T}
    center, axis, radius, length = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8])
    catch
        error("Error parsing 'inside cylinder' constraint data for $filename.")
    end
    return Cylinder{Inside,T}(; center, axis, radius, length)
end

function parse_constraint(::Val{Symbol("outside cylinder")}, data::AbstractVector{<:AbstractString}, ::Type{T}, filename::String) where {T}
    center, axis, radius, length = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8])
    catch
        error("Error parsing 'outside cylinder' constraint data for $filename.")
    end
    return Cylinder{Outside,T}(; center, axis, radius, length)
end
