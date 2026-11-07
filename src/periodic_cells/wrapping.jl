#
# Wrapping of positions and centers of mass into the periodic cell centered at
# `unitcell_center`.
#

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
# The fractional coordinates of a center of mass on a face of the cell (where the
# optimizer projects the molecules pushed against a bound) are computed with a rounding
# error: 0.5000000000000001 for the upper face of a cell of 90.67 Å, for example. The
# test of the closed cell accepts that error, or the molecule would be evaluated on the
# opposite face: in a water slab that ends at a face of the cell, every water molecule
# pushed against that face was wrapped to the region of the other slab, and the
# optimization stalled (no step accepted in a whole chunk).
#
# For dodecahedral and octahedral cells (`wigner_seitz = true`), the image is the one
# closest to the center, in the rhombic dodecahedron or truncated octahedron (the
# Wigner-Seitz cell) where the constraints are defined: a center of mass already in it
# is not moved.
#
function _wrap_cm(cm::SVector{D,T}, unitcell::AbstractMatrix, center::SVector{D,T}, wigner_seitz::Bool) where {D,T}
    wigner_seitz && return _wigner_seitz_image(cm, unitcell, center)
    uc = SMatrix{D,D,T}(unitcell)
    frac = uc \ (cm - center)
    all(f -> abs(f) <= T(0.5) + 8 * eps(T), frac) && return cm
    return wrap_to_center(cm, unitcell, center)
end
_wigner_seitz_image(cm::SVector{3,T}, unitcell::AbstractMatrix, center::SVector{3,T}) where {T} =
    _nearest_periodic_image(cm, unitcell, center)
_wigner_seitz_image(cm::SVector{D,T}, unitcell::AbstractMatrix, center::SVector{D,T}) where {D,T} =
    wrap_to_center(cm, unitcell, center)

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
