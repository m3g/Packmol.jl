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

#
# Build a mapping from molecule index to the index of its first atom (atoms
# of a molecule are contiguous: iat_first:iat_first+natoms-1). This lets
# subset routines (e.g. adjust_constraints!'s active-set loop) address an
# arbitrary, non-contiguous list of molecules in O(1) per molecule instead
# of walking structure-type-contiguous ranges.
#
function _build_mol_iat_first(packmol_system::PackmolSystem)
    mol_iat_first = Vector{Int}(undef, packmol_system.nmols)
    imol = 0
    iat_offset = 0
    for st in packmol_system.structure_types
        for i in 1:st.number_of_molecules
            imol += 1
            mol_iat_first[imol] = iat_offset + (i - 1) * st.natoms + 1
        end
        iat_offset += st.number_of_molecules * st.natoms
    end
    return mol_iat_first
end
