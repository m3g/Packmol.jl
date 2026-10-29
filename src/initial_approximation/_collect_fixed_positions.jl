#
# Collect Cartesian positions of all fixed atoms.
# Returns an empty vector if avoid_overlap is false or there are no fixed atoms.
# The reference coordinates of a fixed molecule already are its atom positions:
# its fixed position and orientation are applied to them when the structure
# is read (`_apply_fixed_center!`), and its molecule position is the identity.
#
function _collect_fixed_positions(packmol_system::PackmolSystem{D,T}) where {D,T}
    if !packmol_system.avoid_overlap
        return SVector{D,T}[]
    end
    fixed_atom_positions = SVector{D,T}[]
    for st in packmol_system.structure_types
        if st.fixed.fixed
            append!(fixed_atom_positions, st.reference_coordinates)
        end
    end
    return fixed_atom_positions
end
