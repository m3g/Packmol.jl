#
# Constraint-only objective function + gradient for pre-optimization.
# Evaluates only constraint penalties (no atom-atom distance penalties), for
# an explicit, possibly non-contiguous list of molecules (`mol_list`), so
# molecules can move freely to satisfy their geometric constraints before the
# main packing optimization.
#
# Used by adjust_constraints!'s active-set loop: molecules already satisfying
# their constraints are geometrically independent of everything else in this
# phase, so excluding them from `mol_list` means this call's cost scales with
# the number of molecules still in violation (plus any just-randomized ones),
# not with the total free-molecule count.
#
function constraint_only_fg_for_mols!(
    g, x,
    fg_output::InteratomicDistanceFG{D,T},
    packmol_system::PackmolSystem{D,T},
    atom_positions::Vector{SVector{D,T}},
    mol_list::Vector{Int},
    mol_structure_type::Vector{Int},
    mol_iat_first::Vector{Int},
) where {D,T}
    x_mol = reinterpret(MoleculePosition{D,T}, x)
    for (k, imol) in enumerate(mol_list)
        packmol_system.molecule_positions[imol] = x_mol[k]
    end
    compute_positions_and_constraints_for_mols!(
        atom_positions, fg_output, packmol_system, mol_list, mol_structure_type, mol_iat_first,
    )
    chain_rule_for_mols!(fg_output, packmol_system, mol_list, mol_structure_type, mol_iat_first)
    g_mol = reinterpret(MoleculePosition{D,T}, g)
    for (k, imol) in enumerate(mol_list)
        g_mol[k] = fg_output.g[imol]
    end
    return fg_output.f
end
