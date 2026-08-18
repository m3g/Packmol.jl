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

# Outside: zero unless the point is simultaneously past the near cap, before
# the far cap, and within the radius (i.e. genuinely inside the forbidden
# cylinder); grows toward the interior. Matches the original Fortran
# Packmol's `outside cylinder` penalty (product of the three violations).
function constraint_penalty(c::Cylinder{Outside}, x)
    _, w, r = _cylinder_axial_radial(c, x)
    s1 = min(-w, zero(w))
    s2 = min(w - c.length, zero(w))
    s3 = min(sum(abs2, r) - c.radius^2, zero(w))
    return c.weight * (s1^2 * s2^2 * s3^2)
end

function constraint_gradient(c::Cylinder{Outside}, x)
    v, w, r = _cylinder_axial_radial(c, x)
    s1 = min(-w, zero(w))
    s2 = min(w - c.length, zero(w))
    s3 = min(sum(abs2, r) - c.radius^2, zero(w))
    return c.weight * (
        -2 * s1 * s2^2 * s3^2 * v +
        2 * s1^2 * s2 * s3^2 * v +
        4 * s1^2 * s2^2 * s3 * r
    )
end

@testitem "Cylinder constructors" begin
    @test InsideCylinder([0,0,0],[0,0,1],1.,2.) == Cylinder{Inside,Float64}([0.,0.,0.],[0.,0.,1.],1.,2.,5.0)
    @test InsideCylinder(center=[0,0,0],axis=[0,0,1],radius=1.,length=2.) == Cylinder{Inside,Float64}([0.,0.,0.],[0.,0.,1.],1.,2.,5.0)
    @test InsideCylinder(center=[0,0,0],axis=[0,0,1],radius=1.,length=2.,weight=2.0) == Cylinder{Inside,Float64}([0.,0.,0.],[0.,0.,1.],1.,2.,2.0)
    @test OutsideCylinder([0,0,0],[0,0,1],1.,2.) == Cylinder{Outside,Float64}([0.,0.,0.],[0.,0.,1.],1.,2.,5.0)
    @test OutsideCylinder(center=[0,0,0],axis=[0,0,1],radius=1.,length=2.) == Cylinder{Outside,Float64}([0.,0.,0.],[0.,0.,1.],1.,2.,5.0)
end

@testitem "Cylinder gradients" begin
    using ForwardDiff
    using StaticArrays
    # Axis-aligned cylinder along z, from z=0 to z=2, radius 1
    c = InsideCylinder(center=[0.0, 0.0, 0.0], axis=[0.0, 0.0, 1.0], radius=1.0, length=2.0)
    for x in (
        SVector(0.3, 0.2, 1.0),   # inside
        SVector(0.3, 0.2, -0.5),  # behind near cap
        SVector(0.3, 0.2, 2.7),   # beyond far cap
        SVector(1.5, 0.2, 1.0),   # outside radius
        SVector(1.5, 0.2, -0.5),  # behind cap and outside radius
    )
        @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    end

    c2 = OutsideCylinder(center=[0.0, 0.0, 0.0], axis=[0.0, 0.0, 1.0], radius=1.0, length=2.0)
    for x in (
        SVector(0.3, 0.2, 1.0),   # inside the forbidden region
        SVector(0.3, 0.2, -0.5),  # behind near cap (already outside)
        SVector(0.3, 0.2, 2.7),   # beyond far cap (already outside)
        SVector(1.5, 0.2, 1.0),   # outside radius (already outside)
    )
        @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c2, x), x) ≈ Packmol.constraint_gradient(c2, x)
    end

    # A tilted, off-center, unnormalized-axis cylinder
    c3 = InsideCylinder(center=[1.0, -2.0, 0.5], axis=[1.0, 1.0, 1.0], radius=0.7, length=3.0)
    for x in (SVector(1.2, -1.7, 1.0), SVector(0.0, 0.0, 0.0), SVector(3.0, 0.0, 3.0))
        @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c3, x), x) ≈ Packmol.constraint_gradient(c3, x)
    end
end

#
# Input parsing functions: must be appended to the "parse_constraint" dictionary:
#
parse_constraint["inside cylinder"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    center, axis, radius, length = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8])
    catch
        error("Error parsing 'inside cylinder' constraint data for $(structure_data[:filename]).")
    end
    return Cylinder{Inside,T}(; center, axis, radius, length)
end

parse_constraint["outside cylinder"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    center, axis, radius, length = try
        parse.(T, data[1:3]), parse.(T, data[4:6]), parse(T, data[7]), parse(T, data[8])
    catch
        error("Error parsing 'outside cylinder' constraint data for $(structure_data[:filename]).")
    end
    return Cylinder{Outside,T}(; center, axis, radius, length)
end
