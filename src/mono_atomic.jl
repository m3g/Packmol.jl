using CellListMap.PeriodicSystems
using StaticArrays
using SPGBox

# Structure carrying the function and gradient of a monoatomic system,
# and the minimum distance between atoms so far.
mutable struct MonoAtomicFG{N,T}
    f::T
    g::Vector{SVector{N,T}}
    dmin::T
end

# Custom copy, reset and reducer functions
import CellListMap.PeriodicSystems: copy_output, reset_output!, reducer
copy_output(x::MonoAtomicFG) = MonoAtomicFG(x.f, copy(x.g), x.dmin)
function reset_output!(output::MonoAtomicFG{N,T}) where {N,T}
    output.f = zero(T)
    fill!(output.g, zero(SVector{N,T}))
    output.dmin = typemax(T)
    return output
end
function reducer(x::MonoAtomicFG, y::MonoAtomicFG)
    x.f += y.f
    x.g .+= y.g
    x.dmin = min(x.dmin, y.dmin)
    return x
end

#=

Updates the function and gradient of the system given a pair of 
particles within the cutoff.

=#
function monoatomic_u_and_g!(x, y, i, j, d2, fg::MonoAtomicFG, tol)
    fg.f += (d2 - tol^2)^2
    g = fg.g
    dv = y - x
    d = sqrt(d2)
    dvdd = 4 * (d2 - tol^2) * dv / d
    g[i] -= dvdd
    g[j] += dvdd
    if d < fg.dmin
        fg.dmin = d
    end
    return fg
end

# Function that computes the function and gradient, and returns the function
# value, to conform with the interface of SPGBox
function fg!(system)
    fg = system.fg
    tol = system.cutoff
    map_pairwise!(
        (x, y, i, j, d2, fg) -> monoatomic_u_and_g!(x, y, i, j, d2, fg, tol),
        system,
    )
    return f
end 

"""
    pack_monoatomic(x::AbstractVector{<:SVector{N,T}}, unitcell, tol)

Pack a monoatomic system with iniital positions `x` and distance tolerance `tol`,
into the unitcell defined by `unitecell`, considering periodic boundary conditions.

The unitcell can be a vector, in the case of Orthorhombic cells, or a matrix, in the
case of triclinic cells.

The coordinates and the unitcells can be two- or three-dimensional. 

"""
function pack_monoatomic(x::AbstractVector{<:SVector{N,T}}, unitcell, tol) where {N,T}
    #  The gradient vector will be allocated by SPGBox, as an auxiliary array
    auxvecs = SPGBox.VAux(x, zero(T))
    # And now we use it as the gradient vector of the periodic system structure
    g = fill!(auxvecs.g, zero(SVector{N,T})) 
    system = PeriodicSystem(
        xpositions = x, 
        unitcell = unitcell,
        cutoff = tol,
        output = MonoAtomicFG(zero(T), g, typemax(T)),
        output_name = :fg,
    )
    spgbox!((g,x) -> fg!(system), x; vaux=auxvecs)
    return x
end

