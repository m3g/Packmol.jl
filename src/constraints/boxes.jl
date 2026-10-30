#
# Box types of constraints: Cubes, Boxes, Planes
#
export Cube, InsideCube, OutsideCube
export Box, InsideBox, OutsideBox

# Default weights
weight_default[:box] = 5.0

#
# Base constraint functions for cubes and boxes. Inside: a wall for each dimension.
# Outside: a point is outside the box if it is outside in any dimension, so the penalty
# is zero unless the point is inside in all of them; then it is the squared distance to
# the nearest face (the shortest way out), and the gradient is along that dimension.
#
function orthogonal_wall(::Type{Inside}, center, side, weight, x)
    xc = x - center
    half = side / 2
    if xc > half
        return weight * (xc - half)^2
    elseif xc < -half
        return weight * (xc + half)^2
    else
        return zero(x)
    end
end

function orthogonal_wall_derivative(::Type{Inside}, center, side, weight, x)
    xc = x - center
    half = side / 2
    if xc > half
        dcdx = 2 * weight * (xc - half)
    elseif xc < -half
        dcdx = 2 * weight * (xc + half)
    else
        dcdx = zero(x)
    end
    return dcdx
end

# Distance to the nearest face of a box (`sides` a vector, or the side of a cube), positive
# for points inside the box in all the dimensions, and the dimension of that face
function _outside_box_depth(center, sides, x)
    xc = x - center
    half(i) = (sides isa Number ? sides : sides[i]) / 2
    depth, k = half(1) - abs(xc[1]), 1
    for i in 2:length(xc)
        d = half(i) - abs(xc[i])
        d < depth && ((depth, k) = (d, i))
    end
    return depth > 0 ? (depth, k) : (zero(depth), 0)
end

function _outside_box_penalty(center, sides, weight, x)
    depth, _ = _outside_box_depth(center, sides, x)
    return weight * depth^2
end

function _outside_box_gradient(center, sides, weight, x)
    depth, k = _outside_box_depth(center, sides, x)
    k == 0 && return zero(x)
    s = x[k] >= center[k] ? one(depth) : -one(depth)
    return typeof(x)(ntuple(i -> i == k ? -2 * weight * depth * s : zero(depth), length(x)))
end

#
# Cube
#
@kwdef struct Cube{Placement,T} <: Constraint{Placement,3,T}
    center::SVector{3,T}
    side::T
    weight::T = weight_default[:box]
end

# Outer constructors that infer T from the given arguments, so that
# `Cube{Placement}(...)` (a partially-applied UnionAll, as used by
# InsideCube/OutsideCube below) can be called without the type
# parameter T needing to be given explicitly. The positional form takes the
# parameters of the input file (`inside cube xmin ymin zmin side`): the lower
# corner and the side; the keyword form, the center and the side.
function Cube{Placement}(corner, side, w=weight_default[:box]; weight=w) where {Placement}
    T = promote_type(eltype(corner), typeof(side), typeof(weight))
    return Cube{Placement}(; center=SVector{3,T}(corner) .+ T(side) / 2, side=T(side), weight=T(weight))
end
function Cube{Placement}(; center, side, weight=weight_default[:box]) where {Placement}
    T = promote_type(eltype(center), typeof(side), typeof(weight))
    return Cube{Placement,T}(SVector{3,T}(center), T(side), T(weight))
end

"""
    InsideCube(lower_corner, side; weight)
    InsideCube(; center, side, weight)

Keeps the atoms inside an axis-aligned cube. With positional arguments, the cube spans
from `lower_corner` to `lower_corner .+ side`, as the input file line
`inside cube xmin ymin zmin side`. With keyword arguments, it is given by its `center`.

# Example

```jldoctest
julia> using Packmol

julia> InsideCube([0.0, 0.0, 0.0], 40.0) == InsideCube(center=[20.0, 20.0, 20.0], side=40.0)
true
```
"""
InsideCube(args...; kargs...) = Cube{Inside}(args...; kargs...)

"""
    OutsideCube(lower_corner, side; weight)
    OutsideCube(; center, side, weight)

Keeps the atoms outside an axis-aligned cube, given as in [`InsideCube`](@ref).
"""
OutsideCube(args...; kargs...) = Cube{Outside}(args...; kargs...)

constraint_penalty(c::Cube{Inside}, x) =
    sum(orthogonal_wall(Inside, c.center[i], c.side, c.weight, x[i]) for i in eachindex(x,c.center))
constraint_gradient(c::Cube{Inside}, x) =
    orthogonal_wall_derivative.(Inside, c.center, c.side, c.weight, x)
constraint_penalty(c::Cube{Outside}, x) = _outside_box_penalty(c.center, c.side, c.weight, x)
constraint_gradient(c::Cube{Outside}, x) = _outside_box_gradient(c.center, c.side, c.weight, x)

#
# Box
#
@kwdef struct Box{Placement,T} <: Constraint{Placement,3,T}
    center::SVector{3, T}
    sides::SVector{3,T}
    weight::T = weight_default[:box]
end

# Outer constructors that infer T from the given arguments, so that
# `Box{Placement}(...)` (a partially-applied UnionAll, as used by
# InsideBox/OutsideBox below) can be called without the type
# parameter T needing to be given explicitly. The positional form takes the
# parameters of the input file (`inside box xmin ymin zmin xmax ymax zmax`):
# the lower and upper corners; the keyword form, the center and the sides.
function Box{Placement}(lower, upper, w=weight_default[:box]; weight=w) where {Placement}
    T = promote_type(eltype(lower), eltype(upper), typeof(weight))
    lo, hi = SVector{3,T}(lower), SVector{3,T}(upper)
    all(lo .<= hi) || throw(ArgumentError("""\n
        The corners of a box must be given as the lower one and the upper one
        (lower <= upper in every dimension).

    """))
    return Box{Placement}(; center=(lo .+ hi) ./ 2, sides=hi .- lo, weight=T(weight))
end
function Box{Placement}(; center, sides, weight=weight_default[:box]) where {Placement}
    T = promote_type(eltype(center), eltype(sides), typeof(weight))
    return Box{Placement,T}(SVector{3,T}(center), SVector{3,T}(sides), T(weight))
end

"""
    InsideBox(lower, upper; weight)
    InsideBox(; center, sides, weight)

Keeps the atoms inside an axis-aligned box. With positional arguments, the box spans
from the `lower` to the `upper` corner, as the input file line
`inside box xmin ymin zmin xmax ymax zmax`. With keyword arguments, it is given by its
`center` and the lengths of its `sides`.

# Example

```jldoctest
julia> using Packmol

julia> InsideBox([0.0, 0.0, 0.0], [40.0, 40.0, 40.0]) == InsideBox(center=[20.0, 20.0, 20.0], sides=[40.0, 40.0, 40.0])
true
```
"""
InsideBox(args...; kargs...) = Box{Inside}(args...; kargs...)

"""
    OutsideBox(lower, upper; weight)
    OutsideBox(; center, sides, weight)

Keeps the atoms outside an axis-aligned box, given as in [`InsideBox`](@ref).
"""
OutsideBox(args...; kargs...) = Box{Outside}(args...; kargs...)

constraint_penalty(c::Box{Inside}, x) =
    sum(orthogonal_wall(Inside, c.center[i], c.sides[i], c.weight, x[i]) for i in eachindex(x,c.center,c.sides))
constraint_gradient(c::Box{Inside}, x) =
    orthogonal_wall_derivative.(Inside, c.center, c.sides, c.weight, x)
constraint_penalty(c::Box{Outside}, x) = _outside_box_penalty(c.center, c.sides, c.weight, x)
constraint_gradient(c::Box{Outside}, x) = _outside_box_gradient(c.center, c.sides, c.weight, x)

#
# Input parsing functions: methods of `parse_constraint`, whose keywords must also be
# listed in `constraint_keywords` (constraint_types.jl):
#
function parse_constraint(::Val{Symbol("inside box")}, data::AbstractVector{<:AbstractString}, ::Type{T}, filename::String) where {T}
    lo, hi = try
        parse.(T, data[1:3]), parse.(T, data[4:6])
    catch
        error("Error parsing 'inside box' constraint data for $filename.")
    end
    center = (lo .+ hi) ./ 2
    sides = hi .- lo
    return Box{Inside,T}(;center, sides)
end

function parse_constraint(::Val{Symbol("outside box")}, data::AbstractVector{<:AbstractString}, ::Type{T}, filename::String) where {T}
    lo, hi = try
        parse.(T, data[1:3]), parse.(T, data[4:6])
    catch
        error("Error parsing 'outside box' constraint data for $filename.")
    end
    center = (lo .+ hi) ./ 2
    sides = hi .- lo
    return Box{Outside,T}(;center, sides)
end

function parse_constraint(::Val{Symbol("inside cube")}, data::AbstractVector{<:AbstractString}, ::Type{T}, filename::String) where {T}
    center, side = try
        parse.(T, data[1:3]), parse.(T, data[4])
    catch
        error("Error parsing 'inside cube' constraint data for $filename.")
    end
    # The data is the lower corner and the side (as in Fortran Packmol)
    return Cube{Inside,T}(; center=center .+ side / 2, side)
end

function parse_constraint(::Val{Symbol("outside cube")}, data::AbstractVector{<:AbstractString}, ::Type{T}, filename::String) where {T}
    center, side = try
        parse.(T, data[1:3]), parse.(T, data[4])
    catch
        error("Error parsing 'outside cube' constraint data for $filename.")
    end
    return Cube{Outside,T}(; center=center .+ side / 2, side)
end
