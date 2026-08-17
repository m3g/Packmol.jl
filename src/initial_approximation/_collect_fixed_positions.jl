#
# Collect Cartesian positions of all fixed atoms.
# Returns an empty vector if avoid_overlap is false or there are no fixed atoms.
#
function _collect_fixed_positions(packmol_system::PackmolSystem{D,T}) where {D,T}
    if !packmol_system.avoid_overlap
        return SVector{D,T}[]
    end
    fixed_atom_positions = SVector{D,T}[]
    for st in packmol_system.structure_types
        if st.fixed.fixed
            mp = st.fixed.position
            R = eulermat(mp.angles)
            for r in st.reference_coordinates
                push!(fixed_atom_positions, R * r + mp.cm)
            end
        end
    end
    return fixed_atom_positions
end
