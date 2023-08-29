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

# Updates the function and gradient of the system given a pair of 
# particles within the cutoff.
function monoatomic_u_and_g!(x, y, i, j, d2, fg::MonoAtomicFG, tol)
    fg.dmin = min(sqrt(d2), fg.dmin)
    fg.f += (d2 - tol^2)^2
    dvdd = 4 * (d2 - tol^2) * (y - x)
    fg.g[i] -= dvdd
    fg.g[j] += dvdd
    return fg
end

# Function that computes the function and gradient, and returns the function
# value and mutates the gradient array, to conform with the interface of SPGBox
function fg!(g, x, system)
    # update positions in the system from x matrix
    for i in eachindex(system.xpositions)
        system.xpositions[i] = eltype(system.xpositions)(@view(x[:, i]))
    end
    tol = system.cutoff
    map_pairwise!(
        (x, y, i, j, d2, fg) -> monoatomic_u_and_g!(x, y, i, j, d2, fg, tol),
        system,
    )
    # update gradient matrix from the gradient vector of vectors
    for i in eachindex(system.fg.g)
        g[:, i] .= system.fg.g[i]
    end
    return system.fg.f
end

function pack_monoatomic_callback(spgresult, system, precision, iprint)
    if spgresult.nit % iprint == 0
        println(
            " Iteration: ", spgresult.nit,
            " Minimum distance: ", system.fg.dmin,
            " Function value: ", spgresult.f
        )
    end
    return false
end

"""
    pack_monoatomic!(positions::AbstractVector{<:SVector{N,T}}, unitcell, tol)

Pack a monoatomic system with iniital positions `x` and distance tolerance `tol`,
into the unitcell defined by `unitecell`, considering periodic boundary conditions.

The unitcell can be a vector, in the case of Orthorhombic cells, or a matrix, in the
case of triclinic cells.

The coordinates and the unitcells can be two- or three-dimensional. 

"""
function pack_monoatomic!(
    positions::AbstractVector{<:SVector{N,T}},
    unitcell,
    tol;
    parallel::Bool=true,
    precision=1e-3,
    iprint=10
) where {N,T}
    #  The gradient vector will be allocated by SPGBox, as an auxiliary array
    x = copy(reinterpret(reshape, T, positions))
    auxvecs = SPGBox.VAux(x, zero(T))
    println("Initializing periodic system...")
    system = PeriodicSystem(
        xpositions=positions,
        unitcell=unitcell,
        cutoff=tol,
        output=MonoAtomicFG(zero(T), similar(positions), typemax(T)),
        output_name=:fg,
        parallel=parallel,
    )
    # spgbox! is called with a single function that computes both the function
    # and the gradient, which in this case is better
    println("Packing...")
    spgboxresult = spgbox!(
        (g, x) -> fg!(g, x, system), x;
        callback=
        (spgresult) -> pack_monoatomic_callback(spgresult, system, precision, iprint),
        vaux=auxvecs,
        nitmax=1000
    )
    # update the positions vector with the new coordinates
    for i in eachindex(positions)
        positions[i] = eltype(positions)(@view(x[:, i]))
    end
    println(spgboxresult)
    println("Minimum distance obtained: ", system.fg.dmin)
    for i in eachindex(positions)
        positions[i] = CellListMap.wrap_to_first(positions[i], system.unitcell)
    end
    return positions
end

@testitem "gradient" begin
    using StaticArrays
    using FiniteDifferences
    using CellListMap.PeriodicSystems
    import Packmol: MonoAtomicFG, fg!

    # Testing function that computes the function value with the definition 
    # of fg! above, to use finite-differences to check the gradient
    function f(x; dimension=2, unitcell=[1, 1], tol=0.1, parallel=false, return_grad=false)
        positions = [SVector{dimension}(x[i:i+dimension-1]) for i in 1:dimension:length(x)]
        system = PeriodicSystem(
            xpositions=positions,
            unitcell=unitcell,
            cutoff=tol,
            output=MonoAtomicFG(0.0, similar(positions), +Inf),
            output_name=:fg,
            parallel=parallel,
        )
        g = similar(x)
        if !return_grad
            return fg!(g, x, system)
        end
        f = fg!(g, x, system)
        return f, g
    end
    x = rand(2, 100)
    @test f(x; return_grad=true)[2] ≈ FiniteDifferences.grad(central_fdm(5, 1), f, x)[1] rtol = 1e-3
    @test f(x; return_grad=true, unitcell=[1 0.5; 0.5 1])[2] ≈
          FiniteDifferences.grad(central_fdm(5, 1), (x) ->
            f(x; unitcell=[1 0.5; 0.5 1]), x)[1] rtol = 1e-3
    x = rand(3, 100)
    @test f(x; dimension=3, unitcell=[1, 1, 1], return_grad=true)[2] ≈
          FiniteDifferences.grad(central_fdm(5, 1), (x) ->
            f(x; dimension=3, unitcell=[1, 1, 1]), x)[1] rtol = 1e-3
    @test f(x; dimension=3, unitcell=[1 0.5 0; 0 1 0.5; 0 0 1], return_grad=true)[2] ≈
          FiniteDifferences.grad(central_fdm(5, 1), (x) ->
            f(x; dimension=3, unitcell=[1 0.5 0; 0 1 0.5; 0 0 1]), x)[1] rtol = 1e-3

end # testitem gradient 



