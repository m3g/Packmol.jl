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
@with_kw struct Cube{Placement,N,T} <: Constraint{Placement,N,T}
    center::SVector{N,T}
    side::T
    weight::T = T(5)
end
function Cube(Placement,center::AbstractVector,side::Number,weight::Number) 
    T = promote_type(eltype(center),typeof(side),typeof(weight)) 
    Cube{Placement,length(center),T}(;
        center=SVector{length(center),T}(center),
        side=T(side),
        weight=T(weight)
    )
end
function Cube(Placement,center::AbstractVector,side::Number) 
    T = promote_type(eltype(center),typeof(side),Int)
    Cube{Placement,length(center),T}(;
        center=SVector{length(center),T}(center),
        side=T(side)
    )
end
InsideCube(args...) = Cube(Inside,args...)
OutsideCube(args...) = Cube(Outside,args...)

constraint_penalty(c::Cube{Placement}, x) where Placement = 
    sum(orthogonal_wall(Placement, c.center[i], c.side, c.weight, x[i]) for i in 1:length(x))
constraint_gradient(c::Cube{Placement}, x) where Placement = 
    orthogonal_wall_derivative.(Placement, c.center, c.side, c.weight, x)

#
# Box
#
@with_kw struct Box{Placement,N,T} <: Constraint{Placement,N,T}
    center::SVector{N,T}
    sides::SVector{N,T}
    weight::T = T(5)
end
function Box(Placement,center::AbstractVector,sides::AbstractVector,weight::Number) 
    @assert length(center) == length(sides) "Box: Center coordinates and sides must have same dimensions."
    N = length(center)
    T = promote_type(eltype(center),eltype(side),typeof(weight)) 
    Box{Placement,N,T}(center=SVector{N,T}(center), sides=SVector{N,T}(sides), weight=T(weight))
end
function Box(Placement,center::AbstractVector,sides::AbstractVector) 
    @assert length(center) == length(sides) "Box: Center coordinates and sides must have same dimensions."
    N = length(center)
    T = promote_type(eltype(center),eltype(sides),Int)
    Box{Placement,N,T}(center=SVector{N,T}(center), sides=SVector{N,T}(sides))
end
InsideBox(args...) = Box(Inside,args...)
OutsideBox(args...) = Box(Outside,args...)

constraint_penalty(c::Box{Placement}, x) where Placement = 
    sum(orthogonal_wall(Placement, c.center[i], c.sides[i], c.weight, x[i]) for i in 1:length(x))
constraint_gradient(c::Box{Placement}, x) where Placement = 
    orthogonal_wall_derivative.(Placement, c.center,c.sides,c.weight,x)

#
# Spheres
#
@with_kw struct Sphere{Placement,N,T} <: Constraint{Placement,N,T}
    center::SVector{N,T}
    radius::T
    weight::T = T(5)
end
function Sphere(Placement,center::AbstractVector,radius::Number,weight::Number) 
    N = length(center) 
    T = promote_type(eltype(center),typeof(radius),typeof(weight)) 
    Sphere{Placement,N,T}(;center=SVector{N,T}(center), radius=T(radius), weight=T(weight))
end
function Sphere(Placement,center::AbstractVector,radius::Number)
    N = length(center) 
    T = promote_type(eltype(center),typeof(radius),Int)
    Sphere{Placement,N,T}(;center=SVector{N,T}(center), radius=T(radius))
end
InsideSphere(args...) = Sphere(Inside, args...)
OutsideSphere(args...) = Sphere(Outside, args...)

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



