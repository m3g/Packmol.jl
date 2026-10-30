using CellListMap: ParticleSystem, NeighborPair, pairwise!
using SPGBox: spgbox!

#
# Structure used to evaluate the function and gradient of the distance between atoms
# in a single pass, using CellListMap.
#
@kwdef mutable struct InteratomicDistanceFG{D,T}
    f::T
    g::Vector{MoleculePosition{D,T}}
    dmin::T
    # Largest overlap `r_i + r_j - d` of a pair, with the radii of the atoms (not the
    # inflated working radii): the convergence criterion
    max_overlap::T
    fmol::Vector{T} # contribution of each molecule to the function value
    # Auxiliary array for gradient: carries the gradient relative to the cartesian coordinates of each atom
    gxcar::Vector{SVector{D,T}}
    # Maximum constraint penalty over all atoms (for convergence checking)
    max_constraint_penalty::T = zero(T)
    # Indices of the atoms (molecules) with a nonzero entry in `gxcar` (`fmol`),
    # set by `cartesian_fg!`. Used by the threaded copies of the output, such that
    # resetting and reducing them costs proportionally to the number of entries
    # actually written, not to the number of atoms (see `reduce_output!` below).
    touched_atoms::Vector{Int32} = Int32[]
    touched_mols::Vector{Int32} = Int32[]
end

function InteratomicDistanceFG{D,T}(packmol_system::PackmolSystem) where {D,T}
    natoms = length(packmol_system.atoms)
    InteratomicDistanceFG(
        f = zero(T),
        g = zero(packmol_system.molecule_positions),
        dmin = typemax(T),
        max_overlap = typemin(T),
        fmol = zeros(T, packmol_system.nmols),
        gxcar = fill(zero(SVector{D,T}), natoms),
        max_constraint_penalty = zero(T),
    )
end

# Custom copy, reset and reducer functions for CellListMap.
#
# The threaded copies (one per batch of the parallel pairwise computation) are
# identified by an empty `g`: the molecule-level gradients are only written by
# chain_rule! after the pairwise computation, so copying them per thread is
# wasted work. The dense `fmol`/`gxcar` arrays of the copies are allocated once,
# zeroed, and kept zeroed except at the entries listed in `touched_mols`/
# `touched_atoms`: resetting a copy, and reducing it into the main output,
# touches only those entries. Near convergence few pairs contribute, while
# zeroing and summing the full arrays of every copy is a memory-bandwidth-bound
# cost proportional to the number of threads times the number of atoms.
_is_thread_copy(x::InteratomicDistanceFG) = isempty(x.g)

function CellListMap.copy_output(x::InteratomicDistanceFG{D,T}) where {D,T}
    InteratomicDistanceFG(
        x.f,
        MoleculePosition{D,T}[],  # empty — not used during pairwise
        x.dmin,
        x.max_overlap,
        zeros(T, length(x.fmol)),
        fill(zero(SVector{D,T}), length(x.gxcar)),
        zero(T),  # max_constraint_penalty — not used by CellListMap threads
        Int32[],
        Int32[],
    )
end
function CellListMap.reset_output!(output::InteratomicDistanceFG{D,T}) where {D,T}
    output.f = zero(T)
    output.dmin = typemax(T)
    output.max_overlap = typemin(T)
    output.max_constraint_penalty = zero(T)
    if _is_thread_copy(output)
        _zero_touched!(output)
    else
        fill!(output.g, zero(MoleculePosition{D,T}))
        fill!(output.fmol, zero(T))
        fill!(output.gxcar, zero(SVector{D,T}))
        empty!(output.touched_atoms)
        empty!(output.touched_mols)
    end
    return output
end
function _zero_touched!(output::InteratomicDistanceFG{D,T}) where {D,T}
    for i in output.touched_atoms
        output.gxcar[i] = zero(SVector{D,T})
    end
    for i in output.touched_mols
        output.fmol[i] = zero(T)
    end
    empty!(output.touched_atoms)
    empty!(output.touched_mols)
    return output
end
function CellListMap.reducer(x::InteratomicDistanceFG, y::InteratomicDistanceFG)
    x.f += y.f
    x.dmin = min(x.dmin, y.dmin)
    x.max_overlap = max(x.max_overlap, y.max_overlap)
    x.fmol .+= y.fmol
    x.gxcar .+= y.gxcar
    return x
end

# The threaded copies are merged into the main output (which was reset before
# the pairwise computation) either sparsely, adding only the touched entries of
# each copy and zeroing them, or densely, summing the full arrays of all copies
# in parallel over index ranges and then zeroing the copies. The sparse merge is
# serial, and is used while the number of touched entries is small relative to
# the cost of the dense merge (number of copies times the number of atoms):
# that is the case near convergence, when few pairs overlap, while early in the
# packing (many overlapping pairs) the dense merge is used.
const _SPARSE_REDUCE_FACTOR = 16
function CellListMap.reduce_output!(
    output::InteratomicDistanceFG{D,T}, output_threaded::Vector{InteratomicDistanceFG{D,T}},
) where {D,T}
    output.f = zero(T)
    output.dmin = typemax(T)
    output.max_overlap = typemin(T)
    ntouched = 0
    for ot in output_threaded
        output.f += ot.f
        output.dmin = min(output.dmin, ot.dmin)
        output.max_overlap = max(output.max_overlap, ot.max_overlap)
        ntouched += length(ot.touched_atoms) + length(ot.touched_mols)
    end
    nmols = length(output.fmol)
    natoms = length(output.gxcar)
    if _SPARSE_REDUCE_FACTOR * ntouched <= length(output_threaded) * (natoms + nmols)
        for ot in output_threaded
            for i in ot.touched_atoms
                output.gxcar[i] += ot.gxcar[i]
            end
            for i in ot.touched_mols
                output.fmol[i] += ot.fmol[i]
            end
            _zero_touched!(ot)
        end
        return output
    end
    @_sync begin
        for (_, mrange) in enumerate(chunks(1:nmols; n=Threads.nthreads()))
            @spawn for i in mrange
                s = zero(T)
                for ot in output_threaded
                    s += ot.fmol[i]
                    ot.fmol[i] = zero(T)
                end
                output.fmol[i] = s
            end
        end
        for (_, arange) in enumerate(chunks(1:natoms; n=Threads.nthreads()))
            @spawn for i in arange
                s = zero(SVector{D,T})
                for ot in output_threaded
                    s += ot.gxcar[i]
                    ot.gxcar[i] = zero(SVector{D,T})
                end
                output.gxcar[i] = s
            end
        end
    end
    for ot in output_threaded
        empty!(ot.touched_atoms)
        empty!(ot.touched_mols)
    end
    return output
end

# Updates the function and gradient of the system given a pair of
# particles within the cutoff.
#
# `atom_radii` implements the same "loose start" heuristic as the original
# Fortran Packmol (getinp.f90's `discale`, default 1.1), generalized to a
# per-atom working radius rather than a single global scale factor: early in
# packing, every atom's own radius (thus its contribution to the pairwise
# target distance) is inflated, giving SPGBox a looser, easier target while
# molecules are far from feasible; each atom's entry in `atom_radii` decays
# toward its own floor (its user-specified or default radius) as improvement
# stalls, and is fattened back up when its molecule is relocated by
# movebad! (see packmol_main.jl). `fg.dmin` tracks the true (unscaled)
# minimum distance and `fg.max_overlap` the largest overlap with the radii of the
# atoms, regardless of `atom_radii`, so convergence checks are unaffected by this
# internal optimization-target adjustment.
function cartesian_fg!(pair::NeighborPair, fg::InteratomicDistanceFG, packmol_system, atom_radii)
    (; x, y, i, j, d2) = pair
    iatom = packmol_system.atoms[i]
    jatom = packmol_system.atoms[j]
    if iatom.molecule_index == jatom.molecule_index
        return fg
    end
    tol = atom_radii[i] + atom_radii[j]
    d = sqrt(d2)
    fg.dmin = min(d, fg.dmin)
    fg.max_overlap = max(fg.max_overlap, iatom.radius + jatom.radius - d)
    if d < tol
        # Function value: add to total function value and to
        # molecular contributions
        fadd = (d - tol)^2
        fg.f += fadd
        _add_fmol!(fg, iatom.molecule_index, fadd)
        _add_fmol!(fg, jatom.molecule_index, fadd)
        # Gradient
        dv = y - x
        dvdd = 2 * (d - tol) * dv / d
        _add_gxcar!(fg, i, -dvdd)
        _add_gxcar!(fg, j, dvdd)
    end
    return fg
end

# Accumulate into `fmol`/`gxcar`, recording the indices of the entries that
# become nonzero (see `touched_atoms`/`touched_mols`). An entry that returns to
# exactly zero and is touched again is listed twice, which is harmless: the
# reduction zeroes each entry of a copy after adding it to the main output.
@inline function _add_fmol!(fg::InteratomicDistanceFG, imol, v)
    fmol = fg.fmol[imol]
    iszero(fmol) && push!(fg.touched_mols, imol)
    fg.fmol[imol] = fmol + v
    return fg
end
@inline function _add_gxcar!(fg::InteratomicDistanceFG, iat, v)
    g = fg.gxcar[iat]
    iszero(g) && push!(fg.touched_atoms, iat)
    fg.gxcar[iat] = g + v
    return fg
end

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
        @_sync for (_, mol_range) in enumerate(chunks(1:nmols_st; n=Threads.nthreads()))
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
    lk = ReentrantLock()
    @_sync for (_, mrange) in enumerate(chunks(1:length(mol_list); n=Threads.nthreads()))
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
                wrapped_cm = HASPBC ? _wrap_cm(cm, unitcell, unitcell_center) : cm
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
            @lock lk begin
                fg_output.f += f_local
                fg_output.max_constraint_penalty = max(fg_output.max_constraint_penalty, max_penalty_local)
            end
        end
    end
    return atom_positions
end

#
# Wrap position to the unit cell image centered at `center`.
# CellListMap.wrap_to_first maps to the cell with one corner at the origin;
# we shift so the cell is centered at `center` instead.
#
# The cell is converted to a static matrix, for the wrapping to be type-stable
# (its dimension is not known from the type of a `Matrix`).
#
function wrap_to_center(x::SVector{D,T}, unitcell::AbstractMatrix, center::SVector{D,T}) where {D,T}
    uc = SMatrix{D,D,T}(unitcell)
    half = SVector{D,T}(ntuple(_ -> T(0.5), D))
    offset = uc * half
    return CellListMap.wrap_to_first(x - center + offset, uc) + center - offset
end

#
# The image of the center of mass of a molecule in which its constraints are
# evaluated (and its atoms are given to the cell lists). A center of mass in the
# closed cell is not moved: the optimizer bounds the centers of mass to the closed
# cell (see `packmol`), and `wrap_to_center` maps the upper faces of the cell to
# the lower ones. A molecule pushed against an upper bound would otherwise be
# evaluated on the opposite face of the cell, where a non-periodic constraint
# (e.g. `over plane`) can be violated by a whole cell length, with a gradient
# that the bound doesn't let it follow: the molecule is trapped.
#
function _wrap_cm(cm::SVector{D,T}, unitcell::AbstractMatrix, center::SVector{D,T}) where {D,T}
    uc = SMatrix{D,D,T}(unitcell)
    frac = uc \ (cm - center)
    all(f -> abs(f) <= T(0.5), frac) && return cm
    return wrap_to_center(cm, unitcell, center)
end

#
# Narrow `packmol_system.unitcell`/`unitcell_center` (typed `Union{Nothing,...}`,
# since PBC is optional) to concretely-typed local values *before* entering a
# @spawn loop that conditionally calls `wrap_to_center`. Reading the Union-typed
# field directly from inside a spawned closure — even guarded by an `if has_pbc`
# check computed outside it — prevents the compiler from proving the wrapped
# position stays a concrete SVector{D,T}, which cascades into heap-boxing every
# constraint_penalty/constraint_gradient call in the loop (measured: ~650x more
# allocation for a mid-size system). The dummy values are never read when
# `has_pbc` is false.
function _typed_unitcell(packmol_system::PackmolSystem{D,T}) where {D,T}
    has_pbc = !isnothing(packmol_system.unitcell)
    unitcell::Matrix{T} = has_pbc ? packmol_system.unitcell : zeros(T, D, D)
    unitcell_center::SVector{D,T} = has_pbc ? packmol_system.unitcell_center : zero(SVector{D,T})
    return has_pbc, unitcell, unitcell_center
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
    @_sync for (_, iatrange) in enumerate(chunks(eachindex(atom_positions); n=Threads.nthreads()))
        @spawn for iat in iatrange
            if HASPBC
                imol = packmol_system.atoms[iat].molecule_index
                cm = packmol_system.molecule_positions[imol].cm
                cm_shift = _wrap_cm(cm, unitcell, unitcell_center) - cm
                cl_xpositions[iat] = atom_positions[iat] + cm_shift
            else
                cl_xpositions[iat] = atom_positions[iat]
            end
        end
    end
    return cl_xpositions
end

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
            @_sync for irange in index_chunks(1:nmols; n=Threads.nthreads())
                @spawn begin
                    f_local = zero(T)
                    max_local = zero(T)
                    for i in irange
                        imol = mol_offset + i
                        R = eulermat(packmol_system.molecule_positions[imol].angles)
                        p = _rotation_fg!(fg_output.gxcar, st, R, at_offset + (i - 1) * st.natoms + 1)
                        fg_output.fmol[imol] += p
                        f_local += p
                        max_local = max(max_local, p)
                    end
                    @lock lk begin
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
    @_sync for iat_range in index_chunks(atoms; n=4 * Threads.nthreads())
        @spawn begin
            fg_local = zero(fg_output.f)
            max_penalty_local = zero(T)
            imol_first = atoms[first(iat_range)].molecule_index
            imol_last = atoms[last(iat_range)].molecule_index
            fmol_first = zero(T)
            fmol_last = zero(T)
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
                    x = x + (_wrap_cm(cm, unitcell, unitcell_center) - cm)
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
            end
            @lock lk_fg begin
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

#
# Combined objective function + gradient for SPGBox.
# Computes pairwise distance penalties + constraint penalties, then
# applies the chain rule to get gradients w.r.t. molecule DOFs.
#
function fg!(g, x,
    cl_system,
    packmol_system::PackmolSystem{D,T},
    atom_positions::Vector{SVector{D,T}},
    free_mol_indices::Vector{Int},
    atom_radii::Vector{T}=map(a -> a.radius, packmol_system.atoms),
) where {D,T}
    # Unpack optimizer variables into free molecule slots
    x_mol = reinterpret(MoleculePosition{D,T}, x)
    @_sync for (_, krange) in enumerate(chunks(1:length(free_mol_indices); n=Threads.nthreads()))
        @spawn for k in krange
            packmol_system.molecule_positions[free_mol_indices[k]] = x_mol[k]
        end
    end
    # Compute Cartesian atomic coordinates from ALL molecule DOFs
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)
    # Update CellListMap positions (wrapped into the canonical cell under
    # PBC — see `_update_cl_positions!` below for why).
    has_pbc, unitcell, unitcell_center = _typed_unitcell(packmol_system)
    if has_pbc
        _update_cl_positions!(cl_system.xpositions, atom_positions, packmol_system, Val(true), unitcell, unitcell_center)
    else
        _update_cl_positions!(cl_system.xpositions, atom_positions, packmol_system, Val(false), unitcell, unitcell_center)
    end
    # Compute pairwise distance penalties and Cartesian gradients
    pairwise!((pair, output) -> cartesian_fg!(pair, output, packmol_system, atom_radii), cl_system,)
    # Add constraint penalties and gradients
    constraint_fg!(cl_system.fg, atom_positions, packmol_system)
    # Zero molecule-level gradients (excluded from CellListMap reset/reduce)
    fill!(cl_system.fg.g, zero(MoleculePosition{D,T}))
    # Chain rule: Cartesian → molecule DOF gradients (for all molecules)
    chain_rule!(cl_system.fg, packmol_system)
    # Pack only free molecule gradients into optimizer gradient vector
    g_mol = reinterpret(MoleculePosition{D,T}, g)
    @_sync for (_, krange) in enumerate(chunks(1:length(free_mol_indices); n=Threads.nthreads()))
        @spawn for k in krange
            g_mol[k] = cl_system.fg.g[free_mol_indices[k]]
        end
    end
    return cl_system.fg.f
end

