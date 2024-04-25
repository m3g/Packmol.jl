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
    if xc > side / 2 || xc < -side / 2
        return weight * (xc - side / 2)^2
    else
        return zero(x)
    end
end

function orthogonal_wall(::Type{Outside}, center, side, weight, x)
    xc = x - center
    if xc < side / 2 && xc > -side / 2
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
struct Cube{Placement,N,T} <: Constraint{Placement,N,T}
    center::SVector{N,T}
    side::T
    weight::T
end
function Cube{Placement}(; center::AbstractVector, side::Number, weight::Number=weight_default[:box]) where {Placement}
    T = promote_type(eltype(center), typeof(side), typeof(weight))
    N = length(center)
    Cube{Placement,N,T}(SVector{N,T}(center), T(side), T(weight))
end
Cube{Placement}(center::AbstractVector, side::Number) where {Placement} =
    Cube{Placement}(; center=center, side=side)
Cube{Placement}(center::AbstractVector, side::Number, weight::Number) where {Placement} =
    Cube{Placement}(; center=center, side=side, weight=weight)

InsideCube(args...; kargs...) = Cube{Inside}(args...; kargs...)
OutsideCube(args...; kargs...) = Cube{Outside}(args...; kargs...)

constraint_penalty(c::Cube{Placement}, x) where {Placement} =
    sum(orthogonal_wall(Placement, c.center[i], c.side, c.weight, x[i]) for i in eachindex(x,c.center))
constraint_gradient(c::Cube{Placement}, x) where {Placement} =
    orthogonal_wall_derivative.(Placement, c.center, c.side, c.weight, x)

#
# Box
#
struct Box{Placement,N,T} <: Constraint{Placement,N,T}
    center::SVector{N,T}
    sides::SVector{N,T}
    weight::T
end
function Box{Placement}(; center::AbstractVector, sides::AbstractVector, weight::Number=weight_default[:box]) where {Placement}
    @assert length(center) == length(sides) "Box: Center coordinates and sides must have same dimensions."
    N = length(center)
    T = promote_type(eltype(center), eltype(sides), typeof(weight))
    Box{Placement,N,T}(SVector{N,T}(center), SVector{N,T}(sides), T(weight))
end
Box{Placement}(center::AbstractVector, sides::AbstractVector) where {Placement} =
    Box{Placement}(; center=center, sides=sides)
Box{Placement}(center::AbstractVector, sides::AbstractVector, weight::Number) where {Placement} =
    Box{Placement}(; center=center, sides=sides, weight=weight)

InsideBox(args...; kargs...) = Box{Inside}(args...; kargs...)
OutsideBox(args...; kargs...) = Box{Outside}(args...; kargs...)

constraint_penalty(c::Box{Placement}, x) where {Placement} =
    sum(orthogonal_wall(Placement, c.center[i], c.sides[i], c.weight, x[i]) for i in eachindex(x,c.center,c.sides))
constraint_gradient(c::Box{Placement}, x) where {Placement} =
    orthogonal_wall_derivative.(Placement, c.center, c.sides, c.weight, x)

#
# Spheres
#
struct Sphere{Placement,N,T} <: Constraint{Placement,N,T}
    center::SVector{N,T}
    radius::T
    weight::T
end
function Sphere{Placement}(; center::AbstractVector, radius::Number, weight::Number=weight_default[:sphere]) where {Placement}
    N = length(center)
    T = promote_type(eltype(center), typeof(radius), typeof(weight))
    Sphere{Placement,N,T}(SVector{N,T}(center), T(radius), T(weight))
end
Sphere{Placement}(center::AbstractVector, radius::Number) where {Placement} =
    Sphere{Placement}(; center=center, radius=radius)

InsideSphere(args...; kargs...) = Sphere{Inside}(args...; kargs...)
OutsideSphere(args...; kargs...) = Sphere{Outside}(args...; kargs...)

function constraint_penalty(c::Sphere{Inside}, x)
    @unpack center, radius, weight = c
    d = norm(x - center)
    if d > radius
        return weight * (d^2 - radius^2)
    else
        return zero(eltype(x))
    end
end

function constraint_penalty(c::Sphere{Outside}, x)
    @unpack center, radius, weight = c
    d = norm(x - center)
    if d < radius
        return weight * (radius^2 - d^2)
    else
        return zero(eltype(x))
    end
end

function constraint_gradient(c::Sphere{Inside}, x)
    @unpack center, radius, weight = c
    dx = x - center
    d = norm(dx)
    if d > radius
        return 2 * weight * dx
    else
        return zero(x)
    end
end

function constraint_gradient(c::Sphere{Outside}, x)
    @unpack center, radius, weight = c
    dx = x - center
    d = norm(dx)
    if d < radius
        return -2 * weight * dx
    else
        return zero(x)
    end
end

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


