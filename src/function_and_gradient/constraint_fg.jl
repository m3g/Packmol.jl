#
# The function and gradient of the constraints: of the positions of the atoms, and of the
# rotations of the molecules (`constrain_rotation`).
#

#
# Add constraint penalties and gradients to the fg output structure.
# When PBC is active, atom positions are wrapped to the unit cell centered
# at unitcell_center before evaluating constraints.
#
function constraint_fg!(
    fg_output::InteratomicDistanceFG{D,T},
    atom_positions::Vector{SVector{D,T}},
    packmol_system::PackmolSystem{D,T},
) where {D,T}
    has_pbc, unitcell, unitcell_center = _typed_unitcell(packmol_system)
    # (`Val(true)`/`Val(false)` in separate branches, rather than `Val(has_pbc)`, for the call
    # to be resolved at compile time.)
    if has_pbc
        _constraint_fg!(fg_output, atom_positions, packmol_system, Val(true), unitcell, unitcell_center)
    else
        _constraint_fg!(fg_output, atom_positions, packmol_system, Val(false), unitcell, unitcell_center)
    end
    return _rotation_constraints_fg!(fg_output, packmol_system)
end

#
# Rotation constraints (`constrain_rotation`, see constraints/rotations.jl) of all
# free molecules: penalties added to `f` and `fmol`, gradients to `gxcar`.
#
function _rotation_constraints_fg!(fg_output::InteratomicDistanceFG{D,T}, packmol_system::PackmolSystem{D,T}) where {D,T}
    lk = ReentrantLock()
    imol_offset = 0
    iat_offset = 0
    for st in packmol_system.structure_types
        nmols = st.fixed.fixed ? 1 : st.number_of_molecules
        if !st.fixed.fixed && !isempty(st.rotation_constraints)
            mol_offset, at_offset = imol_offset, iat_offset
            @_sync for irange in index_chunks(1:nmols; n=_nbatches())
                @spawn begin
                    f_local = zero(T)
                    max_local = zero(T)
                    unconverged_local = Int32[]
                    for i in irange
                        imol = mol_offset + i
                        R = eulermat(packmol_system.molecule_positions[imol].angles)
                        p = _rotation_fg!(fg_output.gxcar, st, R, at_offset + (i - 1) * st.natoms + 1)
                        fg_output.fmol[imol] += p
                        f_local += p
                        max_local = max(max_local, p)
                        if fg_output.count_overlaps[] && p >= fg_output.constraint_precision
                            push!(unconverged_local, imol)
                        end
                    end
                    @lock lk begin
                        append!(fg_output.unconverged_mols, unconverged_local)
                        fg_output.f += f_local
                        fg_output.max_constraint_penalty = max(fg_output.max_constraint_penalty, max_local)
                    end
                end
            end
        end
        imol_offset += nmols
        iat_offset += nmols * st.natoms
    end
    return fg_output
end

# See the comment on _compute_positions_and_constraints_for_mols! above: HASPBC
# must be a type parameter (via Val), not a runtime Bool, so the non-PBC
# specialization never compiles the (allocating) wrap_to_center branch at all.
function _constraint_fg!(
    fg_output::InteratomicDistanceFG{D,T},
    atom_positions::Vector{SVector{D,T}},
    packmol_system::PackmolSystem{D,T},
    ::Val{HASPBC},
    unitcell::Matrix{T},
    unitcell_center::SVector{D,T},
) where {D,T,HASPBC}
    # The atoms are stored molecule by molecule, thus a chunk of atoms spans a
    # contiguous range of molecules, and only the first and last molecules of the
    # range can have atoms in other chunks. The penalties of the molecules
    # strictly inside the range are added directly to `fmol`, and those of the
    # first and last molecules are accumulated locally and added under the lock.
    # More chunks than threads balance the load between threads of different speeds.
    lk_fg = ReentrantLock()
    atoms = packmol_system.atoms
    @_sync for iat_range in index_chunks(atoms; n=4 * _nbatches())
        @spawn begin
            fg_local = zero(fg_output.f)
            max_penalty_local = zero(T)
            imol_first = atoms[first(iat_range)].molecule_index
            imol_last = atoms[last(iat_range)].molecule_index
            fmol_first = zero(T)
            fmol_last = zero(T)
            unconverged_local = Int32[]
            for iat in iat_range
                atom = packmol_system.atoms[iat]
                x = atom_positions[iat]
                if HASPBC
                    # Shift by exactly the offset that wraps this atom's own
                    # *molecule* CM into the canonical cell, rather than
                    # wrapping the atom's own absolute position independently
                    # — see the comment on `_compute_positions_and_constraints_for_mols!`
                    # above for why: an atom offset from a comfortably
                    # mid-cell CM by rotation can independently land right at
                    # a wrap seam even though its own CM never gets near one,
                    # turning a harmless boundary touch into a spurious
                    # jump-to-the-opposite-side-of-the-cell constraint
                    # violation. Shifting every atom of a molecule by the
                    # same amount (its own CM's wrap offset) instead means an
                    # atom only ever crosses a wrap seam when its molecule's
                    # CM does.
                    cm = packmol_system.molecule_positions[atom.molecule_index].cm
                    x = x + (_wrap_cm(cm, unitcell, unitcell_center, packmol_system.wigner_seitz_cell) - cm)
                end
                st = packmol_system.structure_types[atom.structure_type_index]
                atom_penalty = zero(T)
                for ic in atom.constraints
                    c = st.constraints[ic]
                    penalty::T = _constraint_dispatch(constraint_penalty, c, x)::T
                    atom_penalty += penalty
                    fg_output.gxcar[iat]::SVector{D,T} += _constraint_dispatch(constraint_gradient, c, x)::SVector{D,T}
                end
                fg_local += atom_penalty
                imol = atom.molecule_index
                if imol == imol_first
                    fmol_first += atom_penalty
                elseif imol == imol_last
                    fmol_last += atom_penalty
                else
                    fg_output.fmol[imol] += atom_penalty
                end
                max_penalty_local = max(max_penalty_local, atom_penalty)
                if fg_output.count_overlaps[] && atom_penalty >= fg_output.constraint_precision
                    push!(unconverged_local, imol)
                end
            end
            @lock lk_fg begin
                append!(fg_output.unconverged_mols, unconverged_local)
                fg_output.f += fg_local
                fg_output.fmol[imol_first] += fmol_first
                if imol_last != imol_first
                    fg_output.fmol[imol_last] += fmol_last
                end
                fg_output.max_constraint_penalty = max(fg_output.max_constraint_penalty, max_penalty_local)
            end
        end
    end
    return fg_output
end
