#
# Compute the total constraint penalty of a molecule placement, summed over
# all atoms and all constraints applying to each atom.
#
function constraint_penalty_sum(
    mp::MoleculePosition{D,T},
    st::StructureType{D,T},
) where {D,T}
    isempty(st.constraints) && return zero(T)
    R = eulermat(mp.angles)
    fx = zero(T)
    for (j, r) in enumerate(st.reference_coordinates)
        x = R * r + mp.cm
        for ic in st.atom_constraints[j]
            c = st.constraints[ic]
            fx += constraint_penalty(c, x)
        end
    end
    return fx
end
