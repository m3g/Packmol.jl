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
    # Number of pairs of atoms overlapping by more than `overlap_precision` (the
    # convergence criterion), for each pair of structure types (i ≤ j): counted only
    # when `count_overlaps[]` is set (once per loop, for the report). The flag is shared
    # with the threaded copies.
    count_overlaps::Base.RefValue{Bool} = Ref(false)
    overlap_precision::T = zero(T)
    overlaps::Matrix{Int} = zeros(Int, 0, 0)
    # Molecules that do not satisfy the convergence criterion (a pair overlapping by
    # `overlap_precision` or more, or an atom with a constraint penalty of
    # `constraint_precision` or more), with repetitions: recorded only when
    # `count_overlaps[]` is set. Used to classify the molecules as bad, for the active
    # subsystem and for movebad! (see `_mark_unconverged!`).
    constraint_precision::T = zero(T)
    unconverged_mols::Vector{Int32} = Int32[]
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
        overlap_precision = T(packmol_system.tolerance_precision),
        overlaps = zeros(Int, length(packmol_system.structure_types), length(packmol_system.structure_types)),
        constraint_precision = T(packmol_system.constraint_precision),
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
        x.count_overlaps,  # shared flag
        x.overlap_precision,
        zero(x.overlaps),
        x.constraint_precision,
        Int32[],
    )
end
function CellListMap.reset_output!(output::InteratomicDistanceFG{D,T}) where {D,T}
    output.f = zero(T)
    output.dmin = typemax(T)
    output.max_overlap = typemin(T)
    output.max_constraint_penalty = zero(T)
    fill!(output.overlaps, 0)
    empty!(output.unconverged_mols)
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
    x.overlaps .+= y.overlaps
    append!(x.unconverged_mols, y.unconverged_mols)
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
    fill!(output.overlaps, 0)
    empty!(output.unconverged_mols)
    ntouched = 0
    for ot in output_threaded
        output.overlaps .+= ot.overlaps
        append!(output.unconverged_mols, ot.unconverged_mols)
        empty!(ot.unconverged_mols)
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
        for (_, mrange) in enumerate(chunks(1:nmols; n=_nbatches()))
            @spawn for i in mrange
                s = zero(T)
                for ot in output_threaded
                    s += ot.fmol[i]
                    ot.fmol[i] = zero(T)
                end
                output.fmol[i] = s
            end
        end
        for (_, arange) in enumerate(chunks(1:natoms; n=_nbatches()))
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
