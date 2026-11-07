#
# The objective function and its gradient, as given to the solver: the distances between
# atoms and the constraints, with the gradient relative to the variables of the molecules.
#

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
    atom_radii::Vector{T}=map(a -> a.radius, packmol_system.atoms),
) where {D,T}
    # Unpack optimizer variables into free molecule slots
    x_mol = reinterpret(MoleculePosition{D,T}, x)
    @_sync for (_, krange) in enumerate(chunks(1:length(free_mol_indices); n=_nbatches()))
        @spawn for k in krange
            packmol_system.molecule_positions[free_mol_indices[k]] = x_mol[k]
        end
    end
    # Compute Cartesian atomic coordinates from ALL molecule DOFs
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)
    # Update CellListMap positions (wrapped into the canonical cell under
    # PBC — see `_update_cl_positions!` below for why).
    has_pbc, unitcell, unitcell_center = _typed_unitcell(packmol_system)
    if has_pbc
        _update_cl_positions!(cl_system.xpositions, atom_positions, packmol_system, Val(true), unitcell, unitcell_center)
    else
        _update_cl_positions!(cl_system.xpositions, atom_positions, packmol_system, Val(false), unitcell, unitcell_center)
    end
    # Compute pairwise distance penalties and Cartesian gradients
    pairwise!((pair, output) -> cartesian_fg!(pair, output, packmol_system, atom_radii), cl_system,)
    # Add constraint penalties and gradients
    constraint_fg!(cl_system.fg, atom_positions, packmol_system)
    # Zero molecule-level gradients (excluded from CellListMap reset/reduce)
    fill!(cl_system.fg.g, zero(MoleculePosition{D,T}))
    # Chain rule: Cartesian → molecule DOF gradients (for all molecules)
    chain_rule!(cl_system.fg, packmol_system)
    # Pack only free molecule gradients into optimizer gradient vector
    g_mol = reinterpret(MoleculePosition{D,T}, g)
    @_sync for (_, krange) in enumerate(chunks(1:length(free_mol_indices); n=_nbatches()))
        @spawn for k in krange
            g_mol[k] = cl_system.fg.g[free_mol_indices[k]]
        end
    end
    return cl_system.fg.f
end

#
# Marks (in `unconverged`, one entry per molecule) the molecules recorded by the last
# evaluation with `count_overlaps[]` set as not satisfying the convergence criterion.
#
function _mark_unconverged!(unconverged::Vector{Bool}, fg::InteratomicDistanceFG)
    fill!(unconverged, false)
    for imol in fg.unconverged_mols
        unconverged[imol] = true
    end
    return unconverged
end
