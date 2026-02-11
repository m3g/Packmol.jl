using CellListMap: ParticleSystem, NeighborPair, pairwise!
using SPGBox: spgbox!

#
# Structure used to evaluate the function and gradient of the distance between atoms
# in a single pass, using CellListMap.
#
@kwdef mutable struct InteratomicDistanceFG{D,T}
    f::T
    g::Vector{MoleculePosition{D,T}}
    dmin::T
    fmol::Vector{T} # contribution of each molecule to the function value
    # Auxiliary array for gradient: carries the gradient relative to the cartesian coordinates of each atom
    gxcar::Vector{SVector{D,T}}
end
function InteratomicDistanceFG{D,T}(packmol_system::PackmolSystem) where {D,T}
    natoms = length(packmol_system.atoms)
    InteratomicDistanceFG(
        f = zero(T),
        g = zero(packmol_system.molecule_positions),
        dmin = typemax(T),
        fmol = zeros(T, packmol_system.nmols),
        gxcar = fill(zero(SVector{D,T}), natoms),
    )
end

# Custom copy, reset and reducer functions for CellListMap
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
    fill!(output.g, zero(MoleculePosition{D,T}))
    output.dmin = typemax(T)
    fill!(output.fmol, zero(T))
    fill!(output.gxcar, zero(SVector{D,T}))
    return output
end
function CellListMap.reducer(x::InteratomicDistanceFG, y::InteratomicDistanceFG)
    x.f += y.f
    x.g .+= y.g
    x.dmin = min(x.dmin, y.dmin)
    x.fmol .+= y.fmol
    x.gxcar .+= y.gxcar
    return x
end

# Updates the function and gradient of the system given a pair of
# particles within the cutoff.
function cartesian_fg!(pair::NeighborPair, fg::InteratomicDistanceFG, packmol_system)
    (; x, y, i, j, d2) = pair
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

#
# Compute Cartesian atomic positions from rigid-body molecule DOFs.
# atom_position[i] = R(angles) * reference_coord[i] + center_of_mass
#
function compute_atom_positions!(
    atom_positions::Vector{SVector{D,T}},
    molecule_positions,
    packmol_system::PackmolSystem{D,T},
) where {D,T}
    iat = 0
    imol = 0
    for st in packmol_system.structure_types
        ref = st.reference_coordinates
        for _ in 1:st.number_of_molecules
            imol += 1
            mp = molecule_positions[imol]
            R = eulermat(mp.angles...)
            cm = mp.cm
            for j in 1:st.natoms
                iat += 1
                atom_positions[iat] = R * ref[j] + cm
            end
        end
    end
    return atom_positions
end

#
# Wrap position to the unit cell image centered at `center`.
# CellListMap.wrap_to_first maps to the cell with one corner at the origin;
# we shift so the cell is centered at `center` instead.
#
function wrap_to_center(x::SVector{D,T}, unitcell::AbstractMatrix, center::SVector{D,T}) where {D,T}
    half = SVector{D,T}(ntuple(_ -> T(0.5), D))
    offset = SVector{D,T}(unitcell * half)
    return SVector{D,T}(CellListMap.wrap_to_first(x - center + offset, unitcell)) + center - offset
end

#
# Add constraint penalties and gradients to the fg output structure.
# When PBC is active, atom positions are wrapped to the unit cell centered
# at unitcell_center before evaluating constraints.
#
function constraint_fg!(
    fg_output::InteratomicDistanceFG{D,T},
    atom_positions::Vector{SVector{D,T}},
    packmol_system::PackmolSystem{D,T},
) where {D,T}
    has_pbc = !isnothing(packmol_system.unitcell)
    for (iat, atom) in enumerate(packmol_system.atoms)
        x = atom_positions[iat]
        if has_pbc
            x = wrap_to_center(x, packmol_system.unitcell, packmol_system.unitcell_center)
        end
        st = packmol_system.structure_types[atom.structure_type_index]
        for ic in atom.constraints
            c = st.constraints[ic]
            fg_output.f += constraint_penalty(c, x)
            fg_output.gxcar[iat] += constraint_gradient(c, x)
        end
    end
    return fg_output
end

#
# Combined objective function + gradient for SPGBox.
# Computes pairwise distance penalties + constraint penalties, then
# applies the chain rule to get gradients w.r.t. molecule DOFs.
#
function fg!(g, x,
    cl_system,
    packmol_system::PackmolSystem{D,T},
    atom_positions::Vector{SVector{D,T}},
    free_mol_indices::Vector{Int},
) where {D,T}
    # Unpack optimizer variables into free molecule slots
    x_mol = reinterpret(MoleculePosition{D,T}, x)
    for (k, imol) in enumerate(free_mol_indices)
        packmol_system.molecule_positions[imol] = x_mol[k]
    end
    # Compute Cartesian atomic coordinates from ALL molecule DOFs
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)
    # Update CellListMap positions
    cl_system.xpositions .= atom_positions
    # Compute pairwise distance penalties and Cartesian gradients
    pairwise!((pair, output) -> cartesian_fg!(pair, output, packmol_system), cl_system,)
    # Add constraint penalties and gradients
    constraint_fg!(cl_system.fg, atom_positions, packmol_system)
    # Chain rule: Cartesian → molecule DOF gradients (for all molecules)
    chain_rule!(cl_system.fg, packmol_system)
    # Pack only free molecule gradients into optimizer gradient vector
    g_mol = reinterpret(MoleculePosition{D,T}, g)
    for (k, imol) in enumerate(free_mol_indices)
        g_mol[k] = cl_system.fg.g[imol]
    end
    return cl_system.fg.f
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


