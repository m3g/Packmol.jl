using CellListMap: ParticleSystem, map_pairwise, map_pairwise!
using SPGBox: spgbox!

# Structure carrying the function and gradient of a monoatomic system,
# and the minimum distance between atoms so far.
mutable struct MonoAtomicFG{N,T}
    f::T
    g::Vector{SVector{N,T}}
    dmin::T
end

# Custom copy, reset and reducer functions
import CellListMap: copy_output, reset_output!, reducer
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
function monoatomic_u_and_g!(x::T, y::T, i, j, d2, fg::MonoAtomicFG, tol) where {T}
    d = sqrt(d2)
    fg.dmin = min(d, fg.dmin)
    if d < tol
        fg.f += (d - tol)^2
        dv = y - x
        if d > 0
            dvdd = 2 * (d - tol) * dv / d
        else
            vrand = rand(T)
            vrand = (0.1 * tol) * vrand / norm(vrand)
            dvdd = 2 * (d - tol) * vrand
        end
        fg.g[i] -= dvdd
        fg.g[j] += dvdd
    end
    return fg
end

# Function that computes the function and gradient, and returns the function
# value and mutates the gradient array, to conform with the interface of SPGBox
function fg!(g, x, system, tol)
    # update positions in the system from x matrix
    for i in eachindex(system.xpositions)
        system.xpositions[i] = eltype(system.xpositions)(@view(x[:, i]))
    end
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

function pack_monoatomic_callback(spgresult, system, tol, iprint)
    if spgresult.nit % iprint == 0
        println(
            " Iteration: ", spgresult.nit,
            " Minimum distance: ", min(system.fg.dmin, system.cutoff),
            " Function value: ", spgresult.f
        )
    end
    if system.fg.dmin > tol
        return true
    end
    return false
end

"""
    pack_monoatomic!(positions::AbstractVector{<:SVector{N,T}}, unitcell, tol)

Pack a monoatomic system with iniital positions `x` and distance tolerance `tol`,
into the unitcell defined by `unitcell`, considering periodic boundary conditions.

The unitcell can be a vector, in the case of orthorhombic cells, or a matrix, in the
case of triclinic cells.

The coordinates and the unitcells can be two- or three-dimensional. 

# Example

```julia-repl
julia> using Packmol, StaticArrays

julia> coordinates = 100 * rand(SVector{3,Float64}, 10^5);

julia> unitcell = [100.0, 100.0, 100.0]; tolerance = 2.0;

julia> pack_monoatomic!(coordinates, unitcell, tolerance)
```

After packing, the `coordinates` array will have updated positions for the
atoms without overlaps. 

"""
function pack_monoatomic!(
    positions::AbstractVector{<:SVector{N,T}},
    unitcell,
    tol;
    parallel::Bool=true,
    iprint=10
) where {N,T}
    #  The gradient vector will be allocated by SPGBox, as an auxiliary array
    x = copy(reinterpret(reshape, T, positions))
    auxvecs = SPGBox.VAux(x, zero(T))
    println("Initializing periodic system...")
    if unitcell isa AbstractMatrix
        volume = det(unitcell)
    else
        volume = prod(unitcell)
    end
    packing_tol = tol + tol/10 # this avoids slow convergence 
    ncells = min(volume / packing_tol^N, length(positions))
    cutoff = (volume / ncells)^(1/N)
    println("Using cell list cutoff: ", cutoff)
    println("Using packing tolerance: ", packing_tol)
    system = ParticleSystem(
        xpositions=positions,
        unitcell=unitcell,
        cutoff=cutoff,
        output=MonoAtomicFG(zero(T), similar(positions), typemax(T)),
        output_name=:fg,
        parallel=parallel,
    )
    # spgbox! is called with a single function that computes both the function
    # and the gradient, which in this case is better
    println("Packing...")
    spgboxresult = spgbox!(
        (g, x) -> fg!(g, x, system, packing_tol), x;
        callback=
        (spgresult) -> pack_monoatomic_callback(spgresult, system, tol, iprint),
        vaux=auxvecs,
        nitmax=1000
    )
    # update the positions vector with the new coordinates
    for i in eachindex(positions)
        positions[i] = eltype(positions)(@view(x[:, i]))
    end
    println(spgboxresult, "\n")
    println("Minimum distance obtained: ", min(system.fg.dmin, system.cutoff))
    for i in eachindex(positions)
        positions[i] = CellListMap.wrap_to_first(positions[i], system.unitcell)
    end
    return positions
end

@testitem "monoatomic gradient" begin
    using StaticArrays
    using FiniteDifferences
    using CellListMap

    # Testing function that computes the function value with the definition 
    # of fg! above, to use finite-differences to check the gradient
    function f(x; dimension=2, unitcell=[1, 1], tol=0.1, parallel=false, return_grad=false)
        positions = [SVector{dimension}(x[i:i+dimension-1]) for i in 1:dimension:length(x)]
        system = ParticleSystem(
            xpositions=positions,
            unitcell=unitcell,
            cutoff=tol,
            output=Packmol.MonoAtomicFG(0.0, similar(positions), +Inf),
            output_name=:fg,
            parallel=parallel,
        )
        g = similar(x)
        if !return_grad
            return Packmol.fg!(g, x, system, tol)
        end
        f = Packmol.fg!(g, x, system, tol)
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



