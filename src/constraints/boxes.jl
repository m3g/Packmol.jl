#
# Box types of constraints: Cubes, Boxes, Planes
#
export Cube, InsideCube, OutsideCube
export Box, InsideBox, OutsideBox

# Default weights
weight_default[:box] = 5.0

#
# Base constraint functions for cubes and boxes
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

function orthogonal_wall(::Type{Outside}, center, side, weight, x)
    xc = x - center
    if -side / 2 < xc < side / 2
        return weight * (xc - side / 2)^2
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

function orthogonal_wall_derivative(::Type{Outside}, center, side, weight, x)
    xc = x - center
    if xc < side / 2 && xc > -side / 2
        if xc < side / 2
            dcdx = 2 * weight * (xc - side / 2)
        elseif xc > -side / 2
            dcdx = 2 * weight * (side / 2 - xc)
        else
            dcdx = zero(x)
        end
        return dcdx
    else
        return zero(x)
    end
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
# parameter T needing to be given explicitly.
function Cube{Placement}(center, side, weight=weight_default[:box]) where {Placement}
    T = promote_type(eltype(center), typeof(side), typeof(weight))
    return Cube{Placement,T}(SVector{3,T}(center), T(side), T(weight))
end
Cube{Placement}(; center, side, weight=weight_default[:box]) where {Placement} =
    Cube{Placement}(center, side, weight)

InsideCube(args...; kargs...) = Cube{Inside}(args...; kargs...)
OutsideCube(args...; kargs...) = Cube{Outside}(args...; kargs...)

constraint_penalty(c::Cube{Placement}, x) where {Placement} =
    sum(orthogonal_wall(Placement, c.center[i], c.side, c.weight, x[i]) for i in eachindex(x,c.center))
constraint_gradient(c::Cube{Placement}, x) where {Placement} =
    orthogonal_wall_derivative.(Placement, c.center, c.side, c.weight, x)

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
# parameter T needing to be given explicitly.
function Box{Placement}(center, sides, weight=weight_default[:box]) where {Placement}
    T = promote_type(eltype(center), eltype(sides), typeof(weight))
    return Box{Placement,T}(SVector{3,T}(center), SVector{3,T}(sides), T(weight))
end
Box{Placement}(; center, sides, weight=weight_default[:box]) where {Placement} =
    Box{Placement}(center, sides, weight)

InsideBox(args...; kargs...) = Box{Inside}(args...; kargs...)
OutsideBox(args...; kargs...) = Box{Outside}(args...; kargs...)

constraint_penalty(c::Box{Placement}, x) where {Placement} =
    sum(orthogonal_wall(Placement, c.center[i], c.sides[i], c.weight, x[i]) for i in eachindex(x,c.center,c.sides))
constraint_gradient(c::Box{Placement}, x) where {Placement} =
    orthogonal_wall_derivative.(Placement, c.center, c.sides, c.weight, x)

#
# Input parsing functions: must be appended to the "parse_constraint" dictionary:
#
parse_constraint["inside box"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    lo, hi = try
        parse.(T, data[1:3]), parse.(T, data[4:6])
    catch
        error("Error parsing 'inside box' constraint data for $(structure_data[:filename]).")
    end
    center = (lo .+ hi) ./ 2
    sides = hi .- lo
    return Box{Inside,T}(;center, sides)
end

parse_constraint["outside box"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    lo, hi = try
        parse.(T, data[1:3]), parse.(T, data[4:6])
    catch
        error("Error parsing 'outside box' constraint data for $(structure_data[:filename]).")
    end
    center = (lo .+ hi) ./ 2
    sides = hi .- lo
    return Box{Outside,T}(;center, sides)
end

parse_constraint["inside cube"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    center, side = try
        parse.(T, data[1:3]), parse.(T, data[4])
    catch
        error("Error parsing 'inside cube' constraint data for $(structure_data[:filename]).")
    end
    return Cube{Inside,T}(;center, side)
end

parse_constraint["outside cube"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    center, side = try
        parse.(T, data[1:3]), parse.(T, data[4])
    catch
        error("Error parsing 'outside cube' constraint data for $(structure_data[:filename]).")
    end
    return Cube{Outside,T}(;center, side)
end
