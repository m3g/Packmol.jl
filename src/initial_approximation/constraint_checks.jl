#
# The constraints of single molecules: whether they are satisfied, their penalty, and the
# function and gradient of the constraints alone (used in the fits of the initial placement).
#

#
# Check if a molecule placement satisfies all constraints of its structure type.
# Returns true if all atom positions have zero constraint penalty.
#
function satisfies_constraints(
    mp::MoleculePosition{D,T},
    st::StructureType{D,T},
    mol_positions::Vector{SVector{D,T}},
) where {D,T}
    isempty(st.constraints) && isempty(st.rotation_constraints) && return true
    R = eulermat(mp.angles)
    _rotation_penalty(st, R) > zero(T) && return false
    for (j, r) in enumerate(st.reference_coordinates)
        x = R * r + mp.cm
        mol_positions[j] = x
        for ic in st.atom_constraints[j]
            c = st.constraints[ic]
            if _constraint_dispatch(constraint_penalty, c, x) > zero(T)
                return false
            end
        end
    end
    return true
end

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
