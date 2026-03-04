# Copy current free-atom positions from the full atom-position array into the buffer.
function _fill_free_atoms!(
    free_atoms::Vector{SVector{D,T}},
    atom_positions::Vector{SVector{D,T}},
    packmol_system::PackmolSystem{D,T},
) where {D,T}
    ifree = 0
    atom_offset = 0
    for st in packmol_system.structure_types
        if !st.fixed.fixed
            for i in 1:st.number_of_molecules
                iat_first = atom_offset + (i - 1) * st.natoms + 1
                for j in 0:st.natoms-1
                    ifree += 1
                    free_atoms[ifree] = atom_positions[iat_first + j]
                end
            end
        end
        atom_offset += st.number_of_molecules * st.natoms
    end
end
