#
# Build a mapping from molecule index to structure type index.
#
function _build_mol_structure_type(packmol_system::PackmolSystem)
    mol_structure_type = Vector{Int}(undef, packmol_system.nmols)
    imol = 0
    for (ist, st) in enumerate(packmol_system.structure_types)
        for _ in 1:st.number_of_molecules
            imol += 1
            mol_structure_type[imol] = ist
        end
    end
    return mol_structure_type
end
