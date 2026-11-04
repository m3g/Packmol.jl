#
# The bounds of the variables of `n` molecules in the constraint fits of the initial
# approximation and of movebad: under periodic boundary conditions, the centers of mass
# in the box of the cell (`_cm_half_extent`, as in the packing loop), so that they don't
# drift along the directions that no constraint bounds (a molecule above a plane, along
# the plane); without PBC, and for the angles, none (infinite bounds, for the calls to
# have a single type).
#
# They replace the implicit confinement of every atom to the cell, which pushed the
# molecules longer than the tolerance away from the faces of the cell in the fits (the
# lipids of a bilayer crowded in the middle of the cell, 15 of 60 lipids of a leaflet in
# the outer quarter of the square, against 26 for a uniform distribution).
#
function _fit_bounds(packmol_system::PackmolSystem{D,T}, n::Int) where {D,T}
    inf = SVector{D,T}(ntuple(_ -> T(Inf), D))
    lo, hi = if isnothing(packmol_system.unitcell)
        -inf, inf
    else
        h = _cm_half_extent(packmol_system)
        c = packmol_system.unitcell_center::SVector{D,T}
        c - h, c + h
    end
    lower = collect(reinterpret(T, fill(MoleculePosition(lo, -inf), n)))
    upper = collect(reinterpret(T, fill(MoleculePosition(hi, inf), n)))
    return lower, upper
end
