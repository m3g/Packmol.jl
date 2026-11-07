#
# The fixed molecules in the initial placement: their atoms, and the overlaps of the free
# molecules with them.
#

#
# Check if a molecule placed at position `mp` overlaps with any fixed atom.
# Returns true if any atom of the molecule is within `tol` of any fixed atom.
# Uses a ParticleSystem2 where ypositions = fixed atoms (cell list never rebuilt)
# and xpositions is updated to the current molecule atoms before each pairwise! call.
# `mol_positions` is a preallocated buffer for the rotated atom coordinates.
#
# `fixed_lo`/`fixed_hi` are the (undilated) bounding box of the fixed atoms and
# `tol` the overlap tolerance: trial molecules whose own bounding box can't
# possibly come within `tol` of the fixed atoms are rejected immediately,
# without ever touching CellListMap. This matters because the free-molecule
# placement region (where `mp` is sampled from) is typically far larger than
# the fixed structure's own footprint (e.g. a small membrane/capsid embedded
# in a much bigger solvent box) — the vast majority of trials have no chance
# of overlapping and hit this cheap early return.
#
function overlaps_fixed(
    mp::MoleculePosition{D,T},
    ref_coords::Vector{SVector{D,T}},
    fixed_sys::S,
    mol_positions::Vector{SVector{D,T}},
    fixed_lo::SVector{D,T},
    fixed_hi::SVector{D,T},
    tol::T,
) where {D,T,S}
    R = eulermat(mp.angles)
    n = length(ref_coords)
    mol_lo = SVector{D,T}(ntuple(_ -> typemax(T), D))
    mol_hi = SVector{D,T}(ntuple(_ -> typemin(T), D))
    for (j, r) in enumerate(ref_coords)
        x = R * r + mp.cm
        mol_positions[j] = x
        mol_lo = min.(mol_lo, x)
        mol_hi = max.(mol_hi, x)
    end
    if any(mol_hi .< fixed_lo .- tol) || any(mol_lo .> fixed_hi .+ tol)
        return false
    end
    # CellListMap's own cutoff may be larger than `tol` (grown to cap the cell
    # count for large/sparse boxes — see reinitialize_with_bounds.jl), so
    # pairs it finds aren't necessarily real overlaps; filter explicitly.
    CellListMap.update!(fixed_sys; xpositions=@view(mol_positions[1:n]))
    pairwise!(
        (pair, out) -> begin
            pair.d < tol && (out.molecule_badness.value[1] += one(T))
            out
        end,
        fixed_sys,
    )
    return fixed_sys.output.molecule_badness.value[1] > zero(T)
end

# Fallback when there are no fixed atoms (no ParticleSystem built)
function overlaps_fixed(
    mp::MoleculePosition{D,T},
    ref_coords::Vector{SVector{D,T}},
    ::Nothing,
    mol_positions::Vector{SVector{D,T}},
    fixed_lo::SVector{D,T},
    fixed_hi::SVector{D,T},
    tol::T,
) where {D,T}
    return false
end

#
# Build a lightweight CellListMap system (nmols=1 output, one trial molecule's
# atoms queried against the fixed atoms at a time) plus the (undilated)
# bounding box of the fixed atoms, for repeated `overlaps_fixed` queries
# against trial molecule placements. Shared by `reinitialize_with_bounds!`
# and `adjust_constraints!`'s movebad loop — both need the exact same kind of
# system, just used with different (parallel vs. serial) calling patterns.
#
# Returns `(nothing, zero, zero)` when there are no fixed atoms (or
# `avoid_overlap` is false — see `_collect_fixed_positions`): callers pass
# that straight through to `overlaps_fixed`, whose `fixed_sys::Nothing`
# method always returns "no overlap", so no extra branching is needed at the
# call site.
#
#
# `avoid_overlap` (default: the system's) is `true` for the relocations of the packing
# loop, which avoid the fixed atoms whatever the input says (see
# `_collect_fixed_positions`).
#
# Under periodic boundary conditions, the system is built with the periodic cell of the
# system, so that the periodic images of the fixed atoms count too: a fixed structure
# crossing a face of the cell (a protein off the center of the cell, a nanotube, a
# membrane) has atoms near the opposite face, through its images. The bounding box of
# the fixed atoms then does not bound where they are found, and is returned infinite, so
# that `overlaps_fixed` never rejects a trial by it. (The unit cell is a static matrix
# in both cases, diagonal without PBC, for `fixed_sys` to have a single concrete type.)
#
function _build_overlap_check_system(
    packmol_system::PackmolSystem{D,T}, tol::T; avoid_overlap::Bool=packmol_system.avoid_overlap,
) where {D,T}
    max_extent = zero(T)
    for st in packmol_system.structure_types
        st.fixed.fixed && continue
        for r in st.reference_coordinates
            max_extent = max(max_extent, norm(r))
        end
    end
    fixed_atom_positions = _collect_fixed_positions(packmol_system; avoid_overlap)
    isempty(fixed_atom_positions) && return nothing, zero(SVector{D,T}), zero(SVector{D,T})
    lo, hi = compute_bounding_box(fixed_atom_positions)
    pbc_unitcell = packmol_system.unitcell
    fixed_unitcell, volume = if !isnothing(pbc_unitcell)
        lo = SVector{D,T}(ntuple(_ -> typemin(T), D))
        hi = SVector{D,T}(ntuple(_ -> typemax(T), D))
        uc = SMatrix{D,D,T}(pbc_unitcell)
        uc, abs(det(uc))
    else
        margin = max_extent + tol
        sides = T(1.2) * ((hi .+ margin) - (lo .- margin))
        SMatrix{D,D,T}(Diagonal(sides)), prod(sides)
    end
    cutoff = _capped_cutoff(volume, tol, length(fixed_atom_positions), D)
    fixed_sys = build_fixed_particle_system(
        packmol_system; nmols=1, parallel=false, unitcell=fixed_unitcell, cutoff, avoid_overlap,
    )
    return fixed_sys, lo, hi
end

#
# Collect Cartesian positions of all fixed atoms.
# Returns an empty vector if `avoid_overlap` is false or there are no fixed atoms.
# `avoid_overlap` defaults to the system's: the input keyword applies to the initial
# approximation only, and the relocations of the packing loop (movebad) pass `true`,
# since a molecule moved into a fixed structure can take many loops to get out of it.
# The reference coordinates of a fixed molecule already are its atom positions:
# its fixed position and orientation are applied to them when the structure
# is read (`_apply_fixed_center!`), and its molecule position is the identity.
#
function _collect_fixed_positions(
    packmol_system::PackmolSystem{D,T}; avoid_overlap::Bool=packmol_system.avoid_overlap,
) where {D,T}
    if !avoid_overlap
        return SVector{D,T}[]
    end
    fixed_atom_positions = SVector{D,T}[]
    for st in packmol_system.structure_types
        if st.fixed.fixed
            append!(fixed_atom_positions, st.reference_coordinates)
        end
    end
    return fixed_atom_positions
end
