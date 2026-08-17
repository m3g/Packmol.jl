#
# Spherical constraints
#
export Sphere, InsideSphere, OutsideSphere

# Default weights
weight_default[:sphere] = 1e-2

#
# Spheres
#
@kwdef struct Sphere{Placement,T} <: Constraint{Placement,3,T}
    center::SVector{3,T}
    radius::T
    weight::T = weight_default[:sphere]
end
InsideSphere(args...; kargs...) = Sphere{Inside}(args...; kargs...)
OutsideSphere(args...; kargs...) = Sphere{Outside}(args...; kargs...)

function constraint_penalty(c::Sphere{Inside}, x)
    (; center, radius, weight) = c
    w = sum(abs2, x - center) - radius^2
    if w > 0
        return weight * w^2
    else
        return zero(eltype(x))
    end
end

function constraint_penalty(c::Sphere{Outside}, x)
    (; center, radius, weight) = c
    w = sum(abs2, x - center) - radius^2
    if w < 0
        return weight * w^2
    else
        return zero(eltype(x))
    end
end

function constraint_gradient(c::Sphere{Inside}, x)
    (; center, radius, weight) = c
    dx = x - center
    w = sum(abs2, dx) - radius^2
    if w > 0
        return 4 * weight * w * dx
    else
        return zero(x)
    end
end

function constraint_gradient(c::Sphere{Outside}, x)
    (; center, radius, weight) = c
    dx = x - center
    w = sum(abs2, dx) - radius^2
    if w < 0
        return 4 * weight * w * dx
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
    # InsideSphere: point outside the sphere (penalty > 0)
    x = SVector{3,Float64}(1.5, 1.0, 0.)
    c = InsideSphere(center=[0.2, 0., 0.1], radius=0.1)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    # InsideSphere: point inside the sphere (penalty = 0)
    x_in = SVector{3,Float64}(0.21, 0.01, 0.11)
    @test Packmol.constraint_penalty(c, x_in) == 0.0
    @test Packmol.constraint_gradient(c, x_in) == zero(x_in)
    # OutsideSphere: point inside the sphere (penalty > 0)
    c2 = OutsideSphere(center=[0.2, 0., 0.1], radius=1.0)
    x_inside = SVector{3,Float64}(0.3, 0.1, 0.2)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c2, x), x_inside) ≈ Packmol.constraint_gradient(c2, x_inside)
    # OutsideSphere: point outside the sphere (penalty = 0)
    x_out = SVector{3,Float64}(5.0, 5.0, 5.0)
    @test Packmol.constraint_penalty(c2, x_out) == 0.0
    @test Packmol.constraint_gradient(c2, x_out) == zero(x_out)
end

#
# Input parsing functions: must be appended to the "parse_constraint" dictionary:
#
parse_constraint["inside sphere"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    center, radius = try
        parse.(T, data[1:3]), parse.(T, data[4])
    catch
        error("Error parsing 'inside sphere' constraint data for $(structure_data[:filename]).")
    end
    return Sphere{Inside,T}(;center, radius)
end

parse_constraint["outside sphere"] = (structure_data, data::Vector{<:AbstractString}; T=Float64) -> begin
    center, radius = try
        parse.(T, data[1:3]), parse.(T, data[4])
    catch
        error("Error parsing 'outside sphere' constraint data for $(structure_data[:filename]).")
    end
    return Sphere{Outside,T}(;center, radius)
end

