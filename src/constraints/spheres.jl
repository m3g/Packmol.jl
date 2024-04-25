#
# Spherical constraints
#
export Sphere, InsideSphere, OutsideSphere

#
# Generic show function for constraints
#
function Base.show(io::IO, ::MIME"text/plain", c::Constraint)
    println(io, typeof(c))
    fnames = fieldnames(typeof(c))
    for i in 1:length(fnames)-1
        println(io, "    $(fnames[i]) = $(getfield(c,fnames[i]))")
    end
    print(io, "    $(fnames[end]) = $(getfield(c,fnames[end]))")
end

# Default weights
weight_default[:sphere] = 5.0

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

@testitem "Sphere constructors" begin
    @test InsideSphere([0,0,0],1.) == Sphere{Inside,3,Float64}([0.,0.,0.],1.,5.0)
    @test InsideSphere(center=[0,0,0],radius=1.) == Sphere{Inside,3,Float64}([0.,0.,0.],1.,5.0)
    @test InsideSphere(center=[0,0,0],radius=1.,weight=2.0) == Sphere{Inside,3,Float64}([0.,0.,0.],1.,2.0)
    @test OutsideSphere([0,0,0],1.) == Sphere{Outside,3,Float64}([0.,0.,0.],1.,5.0)
    @test OutsideSphere(center=[0,0,0],radius=1.) == Sphere{Outside,3,Float64}([0.,0.,0.],1.,5.0)
    @test OutsideSphere(center=[0,0,0],radius=1.,weight=2.0) == Sphere{Outside,3,Float64}([0.,0.,0.],1.,2.0)
end

@testitem "Sphere gradients" begin
    using ForwardDiff
    using StaticArrays
    x = SVector{3,Float64}(1.5,1.0,0.)
    c = InsideSphere([0.2, 0., 0.1], 0.1)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
    c = OutsideSphere([0.2, 0., 0.1], 1.0)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c,x), x) ≈ Packmol.constraint_gradient(c,x)
end



