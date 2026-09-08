#=
Rhombic dodecahedron periodic box: `pbc dodecahedral cx cy cz d` (input file)
or the `dodecahedral_unitcell` helper (Julia API) build the very same
triclinic unit cell the general `unitcell` keyword does — just computed from
a single "box size" `d` instead of six explicit CRYST1 parameters — so
everything else (CellListMap PBC, implicit confinement, CRYST1 output) needs
no changes at all to support it.

This file also provides `triclinic_to_dodecahedral`/`dodecahedral_to_triclinic`,
a pair of coordinate-remapping helpers between the cell's two equally valid
periodic-image conventions: the ordinary skewed-parallelepiped shape (what
`wrap_to_center` produces, and what CRYST1/PDB output implies) and the
rhombic-dodecahedron-shaped (Wigner-Seitz) "compact" shape, which is more
useful for visualization or analysis centered on a solute (the same
distinction GROMACS draws between `trjconv -ur rect` and `-ur compact`).
=#

export dodecahedral_unitcell
export dodecahedral_to_triclinic, triclinic_to_dodecahedral

"""
    dodecahedral_unitcell(::Type{T}, d) where {T}

3x3 unit cell matrix (columns are the cell vectors) for a rhombic
dodecahedron of "size" `d`, using the standard MD convention (matching
GROMACS's `editconf -bt dodecahedron -d d`):

    v1 = (d, 0, 0)
    v2 = (0, d, 0)
    v3 = (d/2, d/2, d/√2)

equivalently `a = b = c = d`, `α = β = 60°`, `γ = 90°` in CRYST1 terms. `d`
is the guaranteed minimum-image distance in every direction — the point of
the shape is that, for a given cell volume, it gets a bigger minimum-image
distance than a cube does, by trading away the cube's wasted corner volume
relative to a sphere (about 70.7% of a cube's volume for the same `d`).

Use as `PackmolSystem(structure_types; unitcell=dodecahedral_unitcell(Float64, d), unitcell_center=[cx,cy,cz], ...)`
— equivalent to the input file's `pbc dodecahedral cx cy cz d`.
"""
dodecahedral_unitcell(::Type{T}, d) where {T} = Matrix{T}(unitcell_matrix(T, d, d, d, T(60), T(60), T(90)))

#=
    parse_pbc_dodecahedral(::Type{T}, values, D::Int) where {T}

Parses the `cx cy cz d` tail of a `pbc dodecahedral cx cy cz d` input file
line into `(unitcell, unitcell_center)`, as consumed by `read_packmol_input`.
=#
function parse_pbc_dodecahedral(::Type{T}, values, D::Int) where {T}
    D == 3 || throw(ArgumentError("pbc dodecahedral is only supported for D=3, got D=$D"))
    length(values) == 4 || throw(ArgumentError(
        "pbc dodecahedral requires 4 values (cx cy cz d), got $(length(values))"
    ))
    cx, cy, cz, d = (_parse_value(T, "pbc", v) for v in values)
    unitcell = dodecahedral_unitcell(T, d)
    unitcell_center = SVector{D,T}(cx, cy, cz)
    return unitcell, unitcell_center
end

"""
    triclinic_to_dodecahedral(x, unitcell, center)
    triclinic_to_dodecahedral(x, packmol_system::PackmolSystem)

Remaps position(s) `x` to the periodic image closest to `center` — the
rhombic-dodecahedron-shaped (Wigner-Seitz) fundamental domain of the cell,
as opposed to the skewed-parallelepiped one `wrap_to_center` produces.
Searches all 27 combinations of the three lattice vectors' -1/0/+1
multiples: brute force, but the robust way to find the true nearest image
in a skewed cell, since the plain parallelepiped wrap doesn't always land
on it. `x` can be a single `SVector{3}` or a vector of them. `unitcell` and
`center` can also be given together as a `PackmolSystem`.

See also [`dodecahedral_to_triclinic`](@ref) for the reverse mapping.
"""
function triclinic_to_dodecahedral(x::SVector{3,T}, unitcell::AbstractMatrix, center::SVector{3,T}) where {T}
    rel = x - center
    v1 = SVector{3,T}(unitcell[1, 1], unitcell[2, 1], unitcell[3, 1])
    v2 = SVector{3,T}(unitcell[1, 2], unitcell[2, 2], unitcell[3, 2])
    v3 = SVector{3,T}(unitcell[1, 3], unitcell[2, 3], unitcell[3, 3])
    best = rel
    best_dist2 = sum(abs2, rel)
    for n1 in -1:1, n2 in -1:1, n3 in -1:1
        (n1 == 0 && n2 == 0 && n3 == 0) && continue
        cand = rel + n1 * v1 + n2 * v2 + n3 * v3
        dist2 = sum(abs2, cand)
        if dist2 < best_dist2
            best = cand
            best_dist2 = dist2
        end
    end
    return best + center
end

triclinic_to_dodecahedral(xs::AbstractVector{<:SVector{3}}, unitcell::AbstractMatrix, center::SVector{3}) =
    [triclinic_to_dodecahedral(x, unitcell, center) for x in xs]

triclinic_to_dodecahedral(x, packmol_system::PackmolSystem) =
    triclinic_to_dodecahedral(x, packmol_system.unitcell, packmol_system.unitcell_center)

"""
    dodecahedral_to_triclinic(x, unitcell, center)
    dodecahedral_to_triclinic(x, packmol_system::PackmolSystem)

The other direction of [`triclinic_to_dodecahedral`](@ref): remaps
position(s) `x` into the skewed-parallelepiped fundamental domain of the
cell (the ordinary `unitcell` representation, e.g. what CRYST1/PDB output
implies). This is just `wrap_to_center` under a name that pairs with
`triclinic_to_dodecahedral` for this box. `x` can be a single `SVector{3}`
or a vector of them; `unitcell` and `center` can also be given together as
a `PackmolSystem`.
"""
dodecahedral_to_triclinic(x::SVector{3,T}, unitcell::AbstractMatrix, center::SVector{3,T}) where {T} =
    wrap_to_center(x, unitcell, center)

dodecahedral_to_triclinic(xs::AbstractVector{<:SVector{3}}, unitcell::AbstractMatrix, center::SVector{3}) =
    [dodecahedral_to_triclinic(x, unitcell, center) for x in xs]

dodecahedral_to_triclinic(x, packmol_system::PackmolSystem) =
    dodecahedral_to_triclinic(x, packmol_system.unitcell, packmol_system.unitcell_center)
