#
# Types of constraints 
#

abstract type Constraint{Placement,N,T} end

struct Inside end
struct Outside end
struct Over end
struct Below end
export Inside, Outside, Over, Below

export Cube, InsideCube, OutsideCube
export Box, InsideBox, OutsideBox
export Sphere, InsideSphere, OutsideSphere

#
# Generic show function for constraints
#
function Base.show(io::IO,::MIME"text/plain",c::Constraint)
    println(io,typeof(c))
    fnames = fieldnames(typeof(c))
    for i in 1:length(fnames)-1
        println(io,"    $(fnames[i]) = $(getfield(c,fnames[i]))")
    end
    print(io,"    $(fnames[end]) = $(getfield(c,fnames[end]))")
end

# Default weights

const weight_default = (
    box = 5.0, 
    sphere = 5.0
)

#
# Base constraint functions for cubes and boxes
#
function orthogonal_wall(::Type{Inside}, center, side, weight, x)
    xc = x - center
    if xc > side/2 || xc < -side/2
        return weight*(xc - side/2)^2
    else
        return zero(x)
    end
end

function orthogonal_wall(::Type{Outside}, center, side, weight, x)
    xc = x - center
    if xc < side/2 && xc > -side/2
        return weight*(xc - side/2)^2
    else
        return zero(x)
    end
end

function orthogonal_wall_derivative(::Type{Inside}, center, side, weight, x)
    xc = x - center 
    if xc > side/2 
        dcdx = 2*weight*(xc-side/2)
    elseif xc < -side/2
        dcdx = 2*weight*(side/2 - xc)
    else
        dcdx = zero(x)
    end
    return dcdx
end

function orthogonal_wall_derivative(::Type{Outside}, center, side, weight, x)
    xc = x - center 
    if xc < side/2 && xc > -side/2
        if xc < side/2 
            dcdx = 2*weight*(xc-side/2)
        elseif xc > -side/2
            dcdx = 2*weight*(side/2 - xc)
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
function Cube{Placement}(;center::AbstractVector, side::Number, weight::Number=weight_default[:box]) where Placement
    T = promote_type(eltype(center),typeof(side),typeof(weight)) 
    N = length(center)
    Cube{Placement,N,T}(SVector{N,T}(center), T(side), T(weight))
end
Cube{Placement}(center::AbstractVector,side::Number) where Placement  =
    Cube{Placement}(;center=center, side=side)
Cube{Placement}(center::AbstractVector,side::Number,weight::Number) where Placement =
    Cube{Placement}(;center=center, side=side, weight=weight)

InsideCube(args...;kargs...) = Cube{Inside}(args...;kargs...)
OutsideCube(args...;kargs...) = Cube{Outside}(args...;kargs...)

constraint_penalty(c::Cube{Placement}, x) where Placement = 
    sum(orthogonal_wall(Placement, c.center[i], c.side, c.weight, x[i]) for i in 1:length(x))
constraint_gradient(c::Cube{Placement}, x) where Placement = 
    orthogonal_wall_derivative.(Placement, c.center, c.side, c.weight, x)

#
# Box
#
struct Box{Placement,N,T} <: Constraint{Placement,N,T}
    center::SVector{N,T}
    sides::SVector{N,T}
    weight::T
end
function Box{Placement}(;center::AbstractVector,sides::AbstractVector,weight::Number=weight_default[:box]) where Placement
    @assert length(center) == length(sides) "Box: Center coordinates and sides must have same dimensions."
    N = length(center)
    T = promote_type(eltype(center),eltype(sides),typeof(weight)) 
    Box{Placement,N,T}(SVector{N,T}(center), SVector{N,T}(sides), T(weight))
end
Box{Placement}(center::AbstractVector,sides::AbstractVector) where Placement = 
    Box{Placement}(;center=center, sides=sides)
Box{Placement}(center::AbstractVector,sides::AbstractVector, weight::Number) where Placement = 
    Box{Placement}(;center=center, sides=sides, weight=weight)

InsideBox(args...;kargs...) = Box{Inside}(args...;kargs...)
OutsideBox(args...;kargs...) = Box{Outside}(args...;kargs...)

constraint_penalty(c::Box{Placement}, x) where Placement = 
    sum(orthogonal_wall(Placement, c.center[i], c.sides[i], c.weight, x[i]) for i in 1:length(x))
constraint_gradient(c::Box{Placement}, x) where Placement = 
    orthogonal_wall_derivative.(Placement, c.center,c.sides,c.weight,x)

#
# Spheres
#
struct Sphere{Placement,N,T} <: Constraint{Placement,N,T}
    center::SVector{N,T}
    radius::T
    weight::T
end
function Sphere{Placement}(;center::AbstractVector,radius::Number,weight::Number=weight_default[:sphere]) where Placement
    N = length(center) 
    T = promote_type(eltype(center),typeof(radius),typeof(weight)) 
    Sphere{Placement,N,T}(SVector{N,T}(center), T(radius), T(weight))
end
Sphere{Placement}(center::AbstractVector,radius::Number) where Placement =
    Sphere{Placement}(;center=SVector{N,T}(center), radius=T(radius))

InsideSphere(args...;kargs...) = Sphere{Inside}(args...;kargs...)
OutsideSphere(args...;kargs...) = Sphere{Outside}(args...;kargs...)

function constraint_penalty(c::Sphere{Inside}, x)
    @unpack center, radius, weight = c
    d = norm(x - center)
    if d > radius
        return weight*(d^2 - radius^2)
    else
        return zero(eltype(x))
    end
end

function constraint_penalty(c::Sphere{Outside}, x)
    @unpack center, radius, weight = c
    d = norm(x - center)
    if d < radius
        return weight*(radius^2 - d^2)
    else
        return zero(eltype(x))
    end
end

function constraint_gradient(c::Sphere{Inside}, x)
    @unpack center, radius, weight = c
    dx = x - center
    d = norm(dx)
    if d > radius
        return 2*weight*dx
    else
        return zero(x)
    end
end

function constraint_gradient(c::Sphere{Outside}, x)
    @unpack center, radius, weight = c
    dx = x - center
    d = norm(dx)
    if d < radius
        return -2*weight*dx
    else
        return zero(x)
    end
end



