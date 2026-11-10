#
# Custom constraints: a region defined by a function of the position of the atoms,
# f(x) ≤ 0, given in Julia code (there is no input file syntax).
#
# The functions are stored as `FunctionWrapper`s of concrete signatures, so that the
# constraint has a concrete type, and fits in the closed `AnyConstraint` union of the
# packing loop (see constraint_types.jl), whatever function the user gives.
#
export CustomConstraint

# Default weight
weight_default[:custom] = 5.0

struct CustomConstraint{Placement,T} <: Constraint{Placement,3,T}
    f::FunctionWrapper{T,Tuple{SVector{3,T}}}
    gradient::FunctionWrapper{SVector{3,T},Tuple{SVector{3,T}}}
    weight::T
end

"""
    CustomConstraint(f; gradient = nothing, weight = 5.0, T = Float64)

A constraint defined by a function of the position: the atoms are kept in the region where
`f(x) ≤ 0`, for `x` the position of an atom (an `SVector{3}`, in Å). Outside the region,
the penalty is `weight * f(x)^2`, so `f` should grow with the distance to the region (as
`f(x) = norm(x) - R`, for a sphere of radius `R`).
Write `f` in Å, not without units (as `(norm(x)/R)^2 - 1`): for large regions, the
gradient of such a function at small violations is below the stopping tolerance of the
optimizer, and the fits of the initial placement do not converge.

`gradient(x)` is the gradient of `f` at `x`. If it is not given, it is computed by central
finite differences (6 evaluations of `f`). It can be computed exactly by automatic
differentiation, for example with `gradient = x -> ForwardDiff.gradient(f, x)`.

It is a constraint of whole molecules (`constraints`) or of some of their atoms
(`atom_constraints`), as the other constraints. There is no input file syntax for it.

# Example

```julia
using Packmol, LinearAlgebra
# Inside a torus of radii 20 and 6 Å around the z axis
torus = CustomConstraint(x -> hypot(hypot(x[1], x[2]) - 20, x[3]) - 6)
water = structure_type("water.pdb"; number = 500, constraints = [torus])
```
"""
function CustomConstraint(f; gradient = nothing, weight = weight_default[:custom], T::Type = Float64)
    fw = FunctionWrapper{T,Tuple{SVector{3,T}}}(_CustomFunction{T}(f))
    g = isnothing(gradient) ? _FiniteDifferenceGradient{T}(f) : _CustomGradient{T}(gradient)
    gw = FunctionWrapper{SVector{3,T},Tuple{SVector{3,T}}}(g)
    return CustomConstraint{Inside,T}(fw, gw, T(weight))
end

# The functions of the user, converted to the types of the wrappers
struct _CustomFunction{T,F}
    f::F
end
_CustomFunction{T}(f::F) where {T,F} = _CustomFunction{T,F}(f)
(c::_CustomFunction{T})(x::SVector{3,T}) where {T} = T(c.f(x))

struct _CustomGradient{T,G}
    gradient::G
end
_CustomGradient{T}(g::G) where {T,G} = _CustomGradient{T,G}(g)
(c::_CustomGradient{T})(x::SVector{3,T}) where {T} = SVector{3,T}(c.gradient(x))

# Central finite differences, with a step relative to the magnitude of the coordinates
struct _FiniteDifferenceGradient{T,F}
    f::F
end
_FiniteDifferenceGradient{T}(f::F) where {T,F} = _FiniteDifferenceGradient{T,F}(f)
function (c::_FiniteDifferenceGradient{T})(x::SVector{3,T}) where {T}
    h = cbrt(eps(T)) * max(one(T), maximum(abs, x))
    return SVector{3,T}(ntuple(Val(3)) do i
        e = SVector{3,T}(ntuple(j -> j == i ? h : zero(T), Val(3)))
        (T(c.f(x + e)) - T(c.f(x - e))) / (2h)
    end)
end

function constraint_penalty(c::CustomConstraint{Inside,T}, x::SVector{3,T}) where {T}
    v = c.f(x)
    return v > zero(T) ? c.weight * v^2 : zero(T)
end

function constraint_gradient(c::CustomConstraint{Inside,T}, x::SVector{3,T}) where {T}
    v = c.f(x)
    return v > zero(T) ? 2 * c.weight * v * c.gradient(x) : zero(x)
end

function Base.show(io::IO, ::MIME"text/plain", c::CustomConstraint{Inside,T}) where {T}
    print(io, "CustomConstraint{$T}: the region f(x) ≤ 0 of a user function, weight = $(c.weight)")
end
