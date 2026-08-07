#
# Check if a molecule placed at position `mp` overlaps with any fixed atom.
# Returns true if any atom of the molecule is within `tol` of any fixed atom.
# Uses a ParticleSystem2 where ypositions = fixed atoms (cell list never rebuilt)
# and xpositions is updated to the current molecule atoms before each pairwise! call.
# `mol_positions` is a preallocated buffer for the rotated atom coordinates.
#
function overlaps_fixed(
    mp::MoleculePosition{D,T},
    ref_coords::Vector{SVector{D,T}},
    fixed_sys::S,
    mol_positions::Vector{SVector{D,T}},
) where {D,T,S}
    R = eulermat(mp.angles)
    n = length(ref_coords)
    for (j, r) in enumerate(ref_coords)
        mol_positions[j] = R * r + mp.cm
    end
    CellListMap.update!(fixed_sys; xpositions=@view(mol_positions[1:n]))
    pairwise!(
        (pair, out) -> begin
            out.molecule_badness.value[1] += one(T)
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
) where {D,T}
    return false
end