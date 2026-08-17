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
