#
# Plane constraints: above plane / below plane
#
# A plane is defined by a normal vector n = (a, b, c) and a scalar d,
# such that the plane equation is: a*x + b*y + c*z = d
#
# "above plane a b c d" means the atom must satisfy n·x >= d  (Over)
# "below plane a b c d" means the atom must satisfy n·x <= d  (Below)
#
export Plane, AbovePlane, BelowPlane

# Default weights
weight_default[:plane] = 5.0

@kwdef struct Plane{Placement,T} <: Constraint{Placement,3,T}
    normal::SVector{3,T}
    d::T
    weight::T = weight_default[:plane]
end

# Outer constructors that infer T from the given arguments, so that
# `Plane{Placement}(...)` (a partially-applied UnionAll, as used by
# AbovePlane/BelowPlane below) can be called without the type
# parameter T needing to be given explicitly.
function Plane{Placement}(normal, d, weight=weight_default[:plane]) where {Placement}
    T = promote_type(eltype(normal), typeof(d), typeof(weight))
    return Plane{Placement,T}(SVector{3,T}(normal), T(d), T(weight))
end
Plane{Placement}(; normal, d, weight=weight_default[:plane]) where {Placement} =
    Plane{Placement}(normal, d, weight)

AbovePlane(args...; kargs...) = Plane{Over}(args...; kargs...)
BelowPlane(args...; kargs...) = Plane{Below}(args...; kargs...)

# Above plane: atom must satisfy n·x >= d, penalty when n·x < d
_op_plane(::Plane{Over}, x) = x < zero(x)

# Below plane: atom must satisfy n·x <= d, penalty when n·x > d
_op_plane(::Plane{Below}, x) = x > zero(x)

# Constraint penalty (dispatches on Over/Below via _op_plane)
function constraint_penalty(c::Plane, x)
    (; normal, d, weight) = c
    v = dot(normal, x) - d
    if _op_plane(c, v)
        return weight * v^2
    else
        return zero(eltype(x))
    end
end

# Constraint gradient (dispatches on Over/Below via _op_plane)
function constraint_gradient(c::Plane, x)
    (; normal, d, weight) = c
    v = dot(normal, x) - d
    if _op_plane(c, v)
        return 2 * weight * v * normal
    else
        return zero(x)
    end
end

@testitem "Plane constructors" begin
    @test AbovePlane([0, 0, 1], 5.0) == Plane{Over,Float64}([0.0, 0.0, 1.0], 5.0, 5.0)
    @test AbovePlane(normal = [0, 0, 1], d = 5.0) == Plane{Over,Float64}([0.0, 0.0, 1.0], 5.0, 5.0)
    @test AbovePlane(normal = [0, 0, 1], d = 5.0, weight = 2.0) == Plane{Over,Float64}([0.0, 0.0, 1.0], 5.0, 2.0)
    @test BelowPlane([0, 0, 1], 5.0) == Plane{Below,Float64}([0.0, 0.0, 1.0], 5.0, 5.0)
    @test BelowPlane(normal = [0, 0, 1], d = 5.0) == Plane{Below,Float64}([0.0, 0.0, 1.0], 5.0, 5.0)
    @test BelowPlane(normal = [0, 0, 1], d = 5.0, weight = 2.0) == Plane{Below,Float64}([0.0, 0.0, 1.0], 5.0, 2.0)
end

@testitem "Plane gradients" begin
    using ForwardDiff
    using StaticArrays
    # Test above plane: n = (0,0,1), d = 5 -> atom must have z >= 5
    # Point below the plane (should have penalty)
    x = SVector{3,Float64}(1.0, 2.0, 3.0)
    c = AbovePlane([0, 0, 1], 5.0)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    # Point above the plane (no penalty)
    x = SVector{3,Float64}(1.0, 2.0, 7.0)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    # Test below plane: n = (0,0,1), d = 5 -> atom must have z <= 5
    # Point above the plane (should have penalty)
    x = SVector{3,Float64}(1.0, 2.0, 7.0)
    c = BelowPlane([0, 0, 1], 5.0)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    # Point below the plane (no penalty)
    x = SVector{3,Float64}(1.0, 2.0, 3.0)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    # Test with tilted plane normal
    x = SVector{3,Float64}(1.5, 1.0, 0.5)
    c = AbovePlane([1, 1, 1] ./ sqrt(3.0), 2.0)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
    c = BelowPlane([1, 1, 1] ./ sqrt(3.0), 0.1)
    @test ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x) ≈ Packmol.constraint_gradient(c, x)
end

#
# Input parsing functions: must be appended to the "parse_constraint" dictionary:
#
parse_constraint["above plane"] = (structure_data, data::Vector{<:AbstractString}; T = Float64) -> begin
    normal, d = try
        parse.(T, data[1:3]), parse(T, data[4])
    catch
        error("Error parsing 'above plane' constraint data for $(structure_data[:filename]).")
    end
    return Plane{Over,T}(; normal, d)
end
# "over" is the original Fortran Packmol keyword; "above" is a Packmol.jl-only synonym.
parse_constraint["over plane"] = parse_constraint["above plane"]

parse_constraint["below plane"] = (structure_data, data::Vector{<:AbstractString}; T = Float64) -> begin
    normal, d = try
        parse.(T, data[1:3]), parse(T, data[4])
    catch
        error("Error parsing 'below plane' constraint data for $(structure_data[:filename]).")
    end
    return Plane{Below,T}(; normal, d)
end
