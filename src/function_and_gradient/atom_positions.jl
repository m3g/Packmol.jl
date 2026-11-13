#
# The Cartesian positions of the atoms from the variables of the molecules (centers of
# mass and Euler angles), and their copy into the cell lists.
#

#
# Compute Cartesian atomic positions from rigid-body molecule DOFs.
# atom_position[i] = R(angles) * reference_coord[i] + center_of_mass
#
function compute_atom_positions!(
    atom_positions::Vector{SVector{D,T}},
    molecule_positions,
    packmol_system::PackmolSystem{D,T},
) where {D,T}
    imol_offset = 0
    iat_offset = 0
    for st in packmol_system.structure_types
        ref = st.reference_coordinates
        natoms_st = st.natoms
        nmols_st = st.number_of_molecules
        st_imol_offset = imol_offset
        st_iat_offset = iat_offset
        @_sync for (_, mol_range) in enumerate(chunks(1:nmols_st; n=_nbatches()))
            @spawn begin
                for i in mol_range
                    imol = st_imol_offset + i
                    iat_first = st_iat_offset + (i - 1) * natoms_st
                    mp = molecule_positions[imol]
                    R = eulermat(mp.angles)
                    cm = mp.cm
                    for j in 1:natoms_st
                        atom_positions[iat_first + j] = R * ref[j] + cm
                    end
                end
            end
        end
        imol_offset += nmols_st
        iat_offset += nmols_st * natoms_st
    end
    return atom_positions
end

#
# Combined single-pass: compute Cartesian atom positions AND evaluate
# constraint penalties/gradients, for an explicit, possibly non-contiguous
# list of molecules (`mol_list`). Used by adjust_constraints!'s active-set
# loop: only molecules that are still in violation (or were just randomized)
# need their positions/constraints recomputed on a given pass — molecules
# already satisfying their constraints are untouched, since they are
# geometrically independent of everything else in this constraint-only phase.
#
# This does not accumulate into fg_output.f/.fmol/.gxcar via a whole-array
# reset + `+=`: since every atom in mol_list is written by exactly one task
# (mol_list is chunked into disjoint ranges) and constraint penalties/
# gradients depend only on that atom's own position, each entry is fully
# determined within this call — so it is assigned directly rather than
# accumulated, and no reset over molecules outside mol_list is needed. Only
# the shared scalars (fg_output.f, .max_constraint_penalty) need a lock
# across tasks.
function compute_positions_and_constraints_for_mols!(
    atom_positions::Vector{SVector{D,T}},
    fg_output::InteratomicDistanceFG{D,T},
    packmol_system::PackmolSystem{D,T},
    mol_list::Vector{Int},
    mol_structure_type::Vector{Int},
    mol_iat_first::Vector{Int},
) where {D,T}
    has_pbc, unitcell, unitcell_center = _typed_unitcell(packmol_system)
    # (`Val(true)`/`Val(false)` in separate branches, rather than `Val(has_pbc)`, for the call
    # to be resolved at compile time.)
    if has_pbc
        return _compute_positions_and_constraints_for_mols!(
            atom_positions, fg_output, packmol_system, mol_list, mol_structure_type, mol_iat_first,
            Val(true), unitcell, unitcell_center,
        )
    else
        return _compute_positions_and_constraints_for_mols!(
            atom_positions, fg_output, packmol_system, mol_list, mol_structure_type, mol_iat_first,
            Val(false), unitcell, unitcell_center,
        )
    end
end

# `HASPBC` is a type parameter (via Val), not a runtime Bool: `if HASPBC` below is
# resolved at compile time, so the non-PBC specialization never compiles the
# wrap_to_center branch at all. This matters because merely *having* that branch
# in the loop — even when never taken at runtime — makes the compiler infer `x`'s
# type as the join of both branches; since wrap_to_center's return type isn't
# perfectly inferred (it goes through a generic, allocating CellListMap path for
# a plain Matrix unitcell), that instability otherwise poisons every downstream
# constraint_penalty/constraint_gradient call in the loop (measured: ~650x more
# allocation for a mid-size, non-PBC system, on the full-system analog of this
# function — see _constraint_fg! below for the same effect measured directly).
function _compute_positions_and_constraints_for_mols!(
    atom_positions::Vector{SVector{D,T}},
    fg_output::InteratomicDistanceFG{D,T},
    packmol_system::PackmolSystem{D,T},
    mol_list::Vector{Int},
    mol_structure_type::Vector{Int},
    mol_iat_first::Vector{Int},
    ::Val{HASPBC},
    unitcell::Matrix{T},
    unitcell_center::SVector{D,T},
) where {D,T,HASPBC}
    fg_output.f = zero(T)
    fg_output.max_constraint_penalty = zero(T)
    # Partial sums per chunk, added in chunk order: the same result whatever the order
    # in which the tasks finish (see _sum_chunks!)
    mol_chunks = chunks(1:length(mol_list); n=_nbatches())
    f_chunks = zeros(T, length(mol_chunks))
    max_chunks = zeros(T, length(mol_chunks))
    @_sync for (ichunk, mrange) in enumerate(mol_chunks)
        @spawn begin
            f_local = zero(T)
            max_penalty_local = zero(T)
            for k in mrange
                imol = mol_list[k]
                ist = mol_structure_type[imol]
                st = packmol_system.structure_types[ist]
                ref = st.reference_coordinates
                natoms_st = st.natoms
                iat_first = mol_iat_first[imol]
                mp = packmol_system.molecule_positions[imol]
                R = eulermat(mp.angles)
                cm = mp.cm
                # Wrap the molecule's CM once (not each atom's own absolute
                # position independently) — see the comment on
                # `_constraint_fg!` below for why this matters.
                wrapped_cm = HASPBC ? _wrap_cm(cm, unitcell, unitcell_center, packmol_system.wigner_seitz_cell) : cm
                has_constraints = !isempty(st.constraints)
                fmol_local = zero(T)
                for j in 1:natoms_st
                    iat = iat_first + j - 1
                    pos = R * ref[j] + cm
                    atom_positions[iat] = pos
                    gxc = zero(SVector{D,T})
                    if has_constraints
                        x = wrapped_cm + R * ref[j]
                        atom = packmol_system.atoms[iat]
                        for ic in atom.constraints
                            c = st.constraints[ic]
                            penalty::T = _constraint_dispatch(constraint_penalty, c, x)::T
                            f_local += penalty
                            fmol_local += penalty
                            gxc += _constraint_dispatch(constraint_gradient, c, x)::SVector{D,T}
                        end
                    end
                    fg_output.gxcar[iat] = gxc
                end
                if !isempty(st.rotation_constraints)
                    p = _rotation_fg!(fg_output.gxcar, st, R, iat_first)
                    f_local += p
                    fmol_local += p
                end
                fg_output.fmol[imol] = fmol_local
                max_penalty_local = max(max_penalty_local, fmol_local)
            end
            f_chunks[ichunk] = f_local
            max_chunks[ichunk] = max_penalty_local
        end
    end
    _sum_chunks!(fg_output, f_chunks, max_chunks)
    return atom_positions
end

#
# Copy `atom_positions` into CellListMap's own position buffer, wrapping
# into the canonical cell (centered at `unitcell_center`) first when PBC is
# active. `atom_positions` are raw Cartesian coordinates computed straight
# from each molecule's *unbounded* optimization-variable center of mass, so
# over the course of optimization an atom whose true position sits near a
# periodic boundary can end up with a raw coordinate that straddles that
# boundary from the unwrapped side (e.g. z barely negative instead of
# barely below the cell's upper edge). A tiny *physical* displacement
# across that boundary is then a near-full-cell jump in the raw coordinate
# CellListMap sees between two consecutive calls — indistinguishable, from
# CellListMap's side, from the atom genuinely having moved a full cell
# width — which produces spuriously wrong neighbor pairs (and from there a
# wildly incorrect pairwise gradient, confirmed by finite differences at
# one such point: an analytical gradient of exactly 0 against a numerical
# one of order 1e9) for that one atom right at the crossing. Wrapping first
# keeps every atom's coordinate within one canonical copy of the cell on
# every call, so a small physical displacement is always a small raw-
# coordinate displacement too.
#
# Every atom of a molecule is shifted by the *same* offset — the one that
# wraps that molecule's own CM into the canonical cell — rather than each
# atom wrapped independently, for the same reason as `_constraint_fg!`
# below: a molecule's CM can sit close enough to a periodic face that some
# of its own atoms (offset from the CM by rotation) fall on the near side
# of the wrap seam while others fall on the far side. Wrapped
# independently, those atoms end up wrapped by a *different* integer
# multiple of the box than their own CM — landing them in a periodic image
# CellListMap's neighbor search now judges relative to the wrong nearby
# atoms, while the constraint gradient (already CM-consistent) still judges
# the same atom relative to the correct one. That mismatch means the
# distance term and the constraint term can pull a molecule near a boundary
# in ways that don't agree even in sign, which reads externally as the
# molecule being unable to work its way back across the seam.
#
# `Val(HASPBC)` (not a runtime `has_pbc` check) for the same reason as
# `_compute_positions_and_constraints_for_mols!` above: keeps the non-PBC
# specialization from ever compiling the (allocating) wrap_to_center branch
# at all.
function _update_cl_positions!(
    cl_xpositions,
    atom_positions::Vector{SVector{D,T}},
    packmol_system::PackmolSystem{D,T},
    ::Val{HASPBC},
    unitcell::Matrix{T},
    unitcell_center::SVector{D,T},
) where {D,T,HASPBC}
    @_sync for (_, iatrange) in enumerate(chunks(eachindex(atom_positions); n=_nbatches()))
        @spawn for iat in iatrange
            if HASPBC
                imol = packmol_system.atoms[iat].molecule_index
                cm = packmol_system.molecule_positions[imol].cm
                cm_shift = _wrap_cm(cm, unitcell, unitcell_center, packmol_system.wigner_seitz_cell) - cm
                cl_xpositions[iat] = atom_positions[iat] + cm_shift
            else
                cl_xpositions[iat] = atom_positions[iat]
            end
        end
    end
    return cl_xpositions
end
