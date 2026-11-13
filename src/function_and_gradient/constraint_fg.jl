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
    imol_offset = 0
    iat_offset = 0
    for st in packmol_system.structure_types
        nmols = st.fixed.fixed ? 1 : st.number_of_molecules
        if !st.fixed.fixed && !isempty(st.rotation_constraints)
            mol_offset, at_offset = imol_offset, iat_offset
            mol_chunks = index_chunks(1:nmols; n=_nbatches())
            f_chunks = zeros(T, length(mol_chunks))
            max_chunks = zeros(T, length(mol_chunks))
            unconverged_chunks = [Int32[] for _ in 1:length(mol_chunks)]
            @_sync for (ichunk, irange) in enumerate(mol_chunks)
                @spawn begin
                    f_local = zero(T)
                    max_local = zero(T)
                    unconverged_local = unconverged_chunks[ichunk]
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
                    f_chunks[ichunk] = f_local
                    max_chunks[ichunk] = max_local
                end
            end
            _sum_chunks!(fg_output, f_chunks, max_chunks)
            for u in unconverged_chunks
                append!(fg_output.unconverged_mols, u)
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
    # first and last molecules are accumulated locally and added after the tasks, in
    # the order of the chunks, as the sums of the chunks (the same result whatever the
    # order in which the tasks finish). More chunks than threads balance the load
    # between threads of different speeds.
    atoms = packmol_system.atoms
    atom_chunks = index_chunks(atoms; n=4 * _nbatches())
    nchunks = length(atom_chunks)
    f_chunks = zeros(T, nchunks)
    max_chunks = zeros(T, nchunks)
    fmol_first_chunks = zeros(T, nchunks)
    fmol_last_chunks = zeros(T, nchunks)
    unconverged_chunks = [Int32[] for _ in 1:nchunks]
    @_sync for (ichunk, iat_range) in enumerate(atom_chunks)
        @spawn begin
            fg_local = zero(fg_output.f)
            max_penalty_local = zero(T)
            imol_first = atoms[first(iat_range)].molecule_index
            imol_last = atoms[last(iat_range)].molecule_index
            fmol_first = zero(T)
            fmol_last = zero(T)
            unconverged_local = unconverged_chunks[ichunk]
            for iat in iat_range
                atom = packmol_system.atoms[iat]
                x = atom_positions[iat]
                if HASPBC
                    # Every atom of a molecule is shifted by the offset that wraps the
                    # center of mass of the molecule into the cell, not wrapped on its
                    # own: an atom of a molecule well inside the cell can lie at the
                    # boundary of the cell, and wrapped alone it would jump to the
                    # opposite side, with a spurious violation of the constraints. An
                    # atom crosses the boundary only when its molecule does.
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
            f_chunks[ichunk] = fg_local
            max_chunks[ichunk] = max_penalty_local
            fmol_first_chunks[ichunk] = fmol_first
            fmol_last_chunks[ichunk] = fmol_last
        end
    end
    _sum_chunks!(fg_output, f_chunks, max_chunks)
    for (ichunk, iat_range) in enumerate(atom_chunks)
        imol_first = atoms[first(iat_range)].molecule_index
        imol_last = atoms[last(iat_range)].molecule_index
        fg_output.fmol[imol_first] += fmol_first_chunks[ichunk]
        if imol_last != imol_first
            fg_output.fmol[imol_last] += fmol_last_chunks[ichunk]
        end
        append!(fg_output.unconverged_mols, unconverged_chunks[ichunk])
    end
    return fg_output
end

#
# Adds the partial sums of the function of the chunks of a threaded loop to `fg.f`, in
# the order of the chunks, so that the result does not depend on the order in which the
# tasks finish (with a lock, the sums were added as each task finished, and runs with
# the same seed gave different packings), and the maximum penalties to
# `fg.max_constraint_penalty`.
#
function _sum_chunks!(fg, f_chunks::Vector{T}, max_chunks::Vector{T}) where {T}
    for ichunk in eachindex(f_chunks, max_chunks)
        fg.f += f_chunks[ichunk]
        fg.max_constraint_penalty = max(fg.max_constraint_penalty, max_chunks[ichunk])
    end
    return fg
end
