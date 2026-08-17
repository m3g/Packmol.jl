#
# Check if a molecule placement satisfies all constraints of its structure type.
# Returns true if all atom positions have zero constraint penalty.
#
function satisfies_constraints(
    mp::MoleculePosition{D,T},
    st::StructureType{D,T},
    mol_positions::Vector{SVector{D,T}},
) where {D,T}
    isempty(st.constraints) && return true
    R = eulermat(mp.angles)
    for (j, r) in enumerate(st.reference_coordinates)
        x = R * r + mp.cm
        mol_positions[j] = x
        for ic in st.atom_constraints[j]
            c = st.constraints[ic]
            if constraint_penalty(c, x) > zero(T)
                return false
            end
        end
    end
    return true
end