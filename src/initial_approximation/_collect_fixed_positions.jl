#
# Collect Cartesian positions of all fixed atoms.
# Returns an empty vector if `avoid_overlap` is false or there are no fixed atoms.
# `avoid_overlap` defaults to the system's: the input keyword applies to the initial
# approximation only, and the relocations of the packing loop (movebad) pass `true`,
# since a molecule moved into a fixed structure can take many loops to get out of it.
# The reference coordinates of a fixed molecule already are its atom positions:
# its fixed position and orientation are applied to them when the structure
# is read (`_apply_fixed_center!`), and its molecule position is the identity.
#
function _collect_fixed_positions(
    packmol_system::PackmolSystem{D,T}; avoid_overlap::Bool=packmol_system.avoid_overlap,
) where {D,T}
    if !avoid_overlap
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
