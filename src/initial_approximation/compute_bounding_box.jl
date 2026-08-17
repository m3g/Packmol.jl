#
# Bounding box computed from actual atom positions.
# Used after molecule placement to determine the region for CellListMap.
#
function compute_bounding_box(atom_positions::Vector{SVector{D,T}}) where {D,T}
    lo = SVector{D,T}(ntuple(_ -> typemax(T), D))
    hi = SVector{D,T}(ntuple(_ -> typemin(T), D))
    for x in atom_positions
        lo = min.(lo, x)
        hi = max.(hi, x)
    end
    return lo, hi
end