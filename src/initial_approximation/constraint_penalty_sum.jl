#
# Compute the total constraint penalty of a molecule placement, summed over
# all atoms and all constraints applying to each atom.
#
function constraint_penalty_sum(
    mp::MoleculePosition{D,T},
    st::StructureType{D,T},
) where {D,T}
    isempty(st.constraints) && isempty(st.rotation_constraints) && return zero(T)
    R = eulermat(mp.angles)
    fx = _rotation_penalty(st, R)
    for (j, r) in enumerate(st.reference_coordinates)
        x = R * r + mp.cm
        for ic in st.atom_constraints[j]
            c = st.constraints[ic]
            fx += _constraint_dispatch(constraint_penalty, c, x)
        end
    end
    return fx
end

#
# Under PBC, the constraints of a molecule are evaluated at the image of its center of
# mass in the cell (`_wrap_cm`), as in the objective function: the centers of mass are
# drawn anywhere within the bounds, and, in dodecahedral and octahedral cells, part of
# the cell (where the constraints are defined) is only reached through periodic images.
#
function constraint_penalty_sum(
    mp::MoleculePosition{D,T},
    st::StructureType{D,T},
    packmol_system::PackmolSystem{D,T},
) where {D,T}
    unitcell, center = packmol_system.unitcell, packmol_system.unitcell_center
    if !isnothing(unitcell) && !isnothing(center)
        mp = MoleculePosition(_wrap_cm(mp.cm, unitcell, center, packmol_system.wigner_seitz_cell), mp.angles)
    end
    return constraint_penalty_sum(mp, st)
end
