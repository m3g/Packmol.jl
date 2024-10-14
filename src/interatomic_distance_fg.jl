using CellListMap: ParticleSystem, map_pairwise, map_pairwise!
using SPGBox: spgbox!

# Structure used to evaluate the function and gradient of the distance between atoms
# in a single pass, using CellListMap.
@kwdef mutable struct InteratomicDistanceFG{D,T}
    f::T
    g::Vector{MoleculePosition{D,T}}
    dmin::T
    fmol::Vector{T} # contribution of each molecule to the function value
    # Auxiliary array for gradient: carries the gradient relative to the cartesian coordinates of each atom
    gxcar::Vector{SVector{D,T}}
    # Gradient, expressed as Vector{MoleculePosition{D,T}}, which can be reinterpted as a Vector{T}
end
function InteratomicDistanceFG{D,T}(packmol_system::PackmolSystem) where {D,T}
    fg = InteratomicDistanceFG(
        f = zero(T), # function value
        g = zero(packmol_system.molecule_positions), # gradient
        dmin = typemax(T), # minimum distance overall
        fmol = fill(zero(T), packmol_system.nmols), # contribution of each molecule to the function value
        gxcar = fill(zero(SVector{D,T}), length(packmol_system.nmols)), # gradient relative to the cartesian coordinates
    )
    return fg
end


# Custom copy, reset and reducer functions
function CellListMap.copy_output(x::InteratomicDistanceFG)
    InteratomicDistanceFG(
        x.f, 
        copy(x.g),
        x.dmin, 
        copy(x.fmol),
        copy(x.gxcar),
    )
end
function CellListMap.reset_output!(output::InteratomicDistanceFG{D,T}) where {D,T}
    output.f = zero(T)
    fill!(output.g, zero(MoleculePositions{D,T})) 
    output.dmin = typemax(T)
    fill!(output.fmol, zero(T))
    fill!(output.gxcar, zero(SVector{D,T})) 
    return output
end
function CellListMap.reducer(x::InteratomicDistanceFG, y::InteratomicDistanceFG)
    x.f += y.f
    x.g += y.g
    x.dmin = max(x.dmin, y.dmin)
    x.fmol .+= y.fmol
    x.gxcar .+= y.gxcar
    return x
end

# Updates the function and gradient of the system given a pair of 
# particles within the cutoff.
function cartesian_fg!(x::T, y::T, i, j, d2, fg::InteratomicDistanceFG, packmol_system) where {T}
    iatom = packmol_system.atoms[i]
    jatom = packmol_system.atoms[j]
    if iatom.molecule_index == jatom.molecule_index
        return fg
    end
    tol = iatom.radius + jatom.radius
    d = sqrt(d2)
    fg.dmin = min(d, fg.dmin)
    if d < tol
        # Function value: add to total function value and to 
        # molecular contributions
        fadd = (d - tol)^2
        fg.f += fadd
        fg.fmol[iatom.molecule_index] += fadd
        fg.fmol[jatom.molecule_index] += fadd
        # Gradient
        dv = y - x
        dvdd = 2 * (d - tol) * dv / d
        fg.gxcar[i] -= dvdd
        fg.gxcar[j] += dvdd
    end
    return fg
end

# Function that computes the function and gradient, and returns the function
# value and mutates the gradient array, to conform with the interface of SPGBox
# This function mutates the system.fg field. 
function fg!(g, x, system::CellListMap.ParticleSystem1, packmol_system::PackmolSystem)
    system.positions .= reinterpret(SVector{D,T}, x)
    # Compute the function value and component of the gradient relative to the cartesian
    # coordinates for each atom
    map_pairwise!(
        (x, y, i, j, d2, output) -> cartesian_fg!(x, y, i, j, d2, output, packmol_system),
        system,
    )
    # Use the chain rule to compute the gradient relative to the rotations
    # and translations of the molecules
    chain_rule!(fg, packmol_system)
    g .= reinterpret(T, fg.g)
    return system.fg.f
end

function spgbox_callback(spgresult, system, tol, iprint)
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
    packmol()

"""
function packmol(
    packmol_system::PackmolSystem{D,T},
    parallel::Bool=true,
    iprint::Int=1,
    nitmax=1000,
) where {D,T}
    positions = packmol_system.molecule_positions
    x = reinterpret(T, positions)
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
        output=zero(T),
        output_name=:fg,
        parallel=parallel,
    )
    # spgbox! is called with a single function that computes both the function
    # and the gradient, which in this case is better
    println("Packing...")
    spgboxresult = spgbox!(
        (g, x) -> fg!(g, x, system, packmol_system), x;
        callback=
        (spgresult) -> spgbox_callback(spgresult, system, tol, iprint),
        vaux=auxvecs,
        nitmax=nitmax,
    )
    # update the positions vector with the new coordinates
    positions .= reinterpret(MoleculePosition{D,T}, x)
    println(spgboxresult, "\n")
    println("Minimum distance obtained: ", min(system.fg.dmin, system.cutoff))
    return # cartesian coordinates
end

@testitem "monoatomic gradient" begin
    using StaticArrays
    using FiniteDifferences
    using CellListMap
    import Packmol: MonoAtomicFG, fg!

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



