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
    if -side / 2 < xc < side / 2
        return zero(x)
    else
        return weight * (xc - side / 2)^2
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
    if xc > side / 2
        dcdx = 2 * weight * (xc - side / 2)
    elseif xc < -side / 2
        dcdx = 2 * weight * (side / 2 - xc)
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
InsideBox(args...; kargs...) = Box{Inside}(args...; kargs...)
OutsideBox(args...; kargs...) = Box{Outside}(args...; kargs...)

constraint_penalty(c::Box{Placement}, x) where {Placement} =
    sum(orthogonal_wall(Placement, c.center[i], c.sides[i], c.weight, x[i]) for i in eachindex(x,c.center,c.sides))
constraint_gradient(c::Box{Placement}, x) where {Placement} =
    orthogonal_wall_derivative.(Placement, c.center, c.sides, c.weight, x)

@testitem "Box constructors" begin
    @test InsideCube([0,0,0],1.) == Cube{Inside,3,Float64}([0.,0.,0.],1.,5.0)
    @test InsideCube(center=[0,0,0],side=1.) == Cube{Inside,3,Float64}([0.,0.,0.],1.,5.0)
    @test InsideCube(center=[0,0,0],side=1.,weight=2.0) == Cube{Inside,3,Float64}([0.,0.,0.],1.,2.0)
    @test OutsideCube([0,0,0],1.) == Cube{Outside,3,Float64}([0.,0.,0.],1.,5.0)
    @test OutsideCube(center=[0,0,0],side=1.) == Cube{Outside,3,Float64}([0.,0.,0.],1.,5.0)
    @test OutsideCube(center=[0,0,0],side=1.,weight=2.0) == Cube{Outside,3,Float64}([0.,0.,0.],1.,2.0)
    @test InsideBox([0,0,0],[1.,1.,1.]) == Box{Inside,3,Float64}([0.,0.,0.],[1.,1.,1.],5.0)
    @test InsideBox(center=[0,0,0],sides=[1.,1.,1.]) == Box{Inside,3,Float64}([0.,0.,0.],[1.,1.,1.],5.0)
    @test InsideBox(center=[0,0,0],sides=[1.,1.,1.],weight=2.0) == Box{Inside,3,Float64}([0.,0.,0.],[1.,1.,1.],2.0)
    @test OutsideBox([0,0,0],[1.,1.,1.]) == Box{Outside,3,Float64}([0.,0.,0.],[1.,1.,1.],5.0)
    @test OutsideBox(center=[0,0,0],sides=[1.,1.,1.]) == Box{Outside,3,Float64}([0.,0.,0.],[1.,1.,1.],5.0)
    @test OutsideBox(center=[0,0,0],sides=[1.,1.,1.],weight=2.0) == Box{Outside,3,Float64}([0.,0.,0.],[1.,1.,1.],2.0)
end

@testitem "Constraint gradients" begin
    using ForwardDiff
    using StaticArrays
    x = SVector{3,Float64}(1.5,1.0,0.)
    c = InsideCube([0.2, 0., 0.1], 0.5)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
    c = OutsideCube([0.2, 0., 0.1], 2.)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
    c = InsideBox([0.2, 0., 0.1], [0.5, 0.7, 1.0])
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
    c = OutsideBox([0.2, 0., 0.1], [0.5, 0.7, 1.0])
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
end

#
# Input parsing functions: must be appended to the "parse_constraint" dictionary:
#
parse_constraint["inside box"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    center, sides = try
        parse.(T, data[1:3]), parse.(T, data[4:6])
    catch
        error("Error parsing 'inside box' constraint data for $(structure_data[:filename]).")
    end
    return Box{Inside,T}(;center, sides)
end

parse_constraint["outside box"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    center, sides = try
        parse.(T, data[1:3]), parse.(T, data[4:6])
    catch
        error("Error parsing 'outside box' constraint data for $(structure_data[:filename]).")
    end
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

