#
# Spherical constraints
#
export Sphere, InsideSphere, OutsideSphere

# Default weights
weight_default[:sphere] = 5.0

#
# Spheres
#
@kwdef struct Sphere{Placement,N,T} <: Constraint{Placement,N,T}
    center::SVector{N,T}
    radius::T
    weight::T = weight_default[:sphere]
end

InsideSphere(args...; kargs...) = Sphere{Inside}(args...; kargs...)
OutsideSphere(args...; kargs...) = Sphere{Outside}(args...; kargs...)

function constraint_penalty(c::Sphere{Inside}, x)
    (; center, radius, weight) = c
    d = norm(x - center)
    if d > radius
        return weight * (d^2 - radius^2)
    else
        return zero(eltype(x))
    end
end

function constraint_penalty(c::Sphere{Outside}, x)
    (; center, radius, weight) = c
    d = norm(x - center)
    if d < radius
        return weight * (radius^2 - d^2)
    else
        return zero(eltype(x))
    end
end

function constraint_gradient(c::Sphere{Inside}, x)
    (; center, radius, weight) = c
    dx = x - center
    d = norm(dx)
    if d > radius
        return 2 * weight * dx
    else
        return zero(x)
    end
end

function constraint_gradient(c::Sphere{Outside}, x)
    (; center, radius, weight) = c
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

#
# Input parsing functions: must be appended to the "parse_constraint" dictionary:
#
parse_constraint["inside sphere"] = (structure_data, data::Vector{<:AbstractString}; T=Float64, N=3) -> begin
    center, radius = try
        parse.(T, data[2+1:2+N]), parse.(T, last(data))
    catch
        error("Error parsing 'inside sphere' constraint data for $(structure_data[:filename]).")
    end
    return Sphere{Inside,N,T}(;center, radius)
end

parse_constraint["outside sphere"] = (structure_data, data::Vector{<:AbstractString}; T=Float64, N=3) -> begin
    center, radius = try
        parse.(T, data[2+1:2+N]), parse.(T, last(data))
    catch
        error("Error parsing 'outside sphere' constraint data for $(structure_data[:filename]).")
    end
    return Sphere{Outside,N,T}(;center, radius)
end

