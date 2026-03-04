#
# Constraint-only objective function + gradient for pre-optimization.
# Evaluates only constraint penalties (no atom-atom distance penalties),
# so molecules can move freely to satisfy their geometric constraints
# before the main packing optimization.
#
function constraint_only_fg!(
    g, x,
    fg_output::InteratomicDistanceFG{D,T},
    packmol_system::PackmolSystem{D,T},
    atom_positions::Vector{SVector{D,T}},
    free_mol_indices::Vector{Int},
) where {D,T}
    # Unpack optimizer variables into free molecule slots
    x_mol = reinterpret(MoleculePosition{D,T}, x)
    for (k, imol) in enumerate(free_mol_indices)
        packmol_system.molecule_positions[imol] = x_mol[k]
    end
    # Reset fg output (no CellListMap to do it for us)
    CellListMap.reset_output!(fg_output)
    # Single pass: compute atom positions AND constraint penalties/gradients
    compute_positions_and_constraints!(atom_positions, fg_output, packmol_system.molecule_positions, packmol_system)
    # Chain rule: Cartesian → molecule DOF gradients (for all molecules)
    chain_rule!(fg_output, packmol_system)
    # Pack only free molecule gradients into optimizer gradient vector
    g_mol = reinterpret(MoleculePosition{D,T}, g)
    for (k, imol) in enumerate(free_mol_indices)
        g_mol[k] = fg_output.g[imol]
    end
    return fg_output.f
end
