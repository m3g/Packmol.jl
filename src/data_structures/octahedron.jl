#=
Truncated octahedron periodic box: `pbc octahedral cx cy cz d` (input file)
or the `octahedral_unitcell` helper (Julia API) build the very same
triclinic unit cell the general `unitcell` keyword does — just computed from
a single "box size" `d` instead of six explicit CRYST1 parameters — so
everything else (CellListMap PBC, implicit confinement, CRYST1 output) needs
no changes at all to support it. Mirrors dodecahedron.jl's own structure —
see that file for the shared `_nearest_periodic_image` "compact" wrap
(periodic_cells.jl) both box types build on.

This file also provides `triclinic_to_octahedral`/`octahedral_to_triclinic`,
a pair of coordinate-remapping helpers between the cell's two equally valid
periodic-image conventions: the ordinary skewed-parallelepiped shape (what
`wrap_to_center` produces, and what CRYST1/PDB output implies) and the
truncated-octahedron-shaped (Wigner-Seitz) "compact" shape, which is more
useful for visualization or analysis centered on a solute (the same
distinction GROMACS draws between `trjconv -ur rect` and `-ur compact`).
Each also has a single-argument `(::PackmolSystem)` method that instead
toggles `packmol_system.periodic_boundary_style` — the flag `get_atoms`
(see write_output.jl) and `write_output` read to decide which shape to wrap
molecule positions into when they materialize atomic coordinates.
=#

export octahedral_unitcell
export octahedral_to_triclinic, triclinic_to_octahedral

"""
    octahedral_unitcell(::Type{T}, d) where {T}

3x3 unit cell matrix (columns are the cell vectors) for a truncated
octahedron of "size" `d`, using the standard MD convention (matching
GROMACS's `editconf -bt octahedron -d d`):

    v1 = (d, 0, 0)
    v2 = (d/3, 2√2 d/3, 0)
    v3 = (-d/3, √2 d/3, √6 d/3)

equivalently `a = b = c = d`, `α = γ = acosd(1/3) ≈ 70.53°`,
`β = acosd(-1/3) ≈ 109.47°` in CRYST1 terms — the Wigner-Seitz cell of a BCC
lattice, built from a Minkowski-reduced (shortest-vectors) basis for it. `d`
is the guaranteed minimum-image distance in every direction. For a given
cell volume this gets a bigger minimum-image distance than a cube does
(about 77.0% of a cube's volume for the same `d`) — a bit more volume than
the rhombic dodecahedron needs for the same `d` (see
[`dodecahedral_unitcell`](@ref), ≈70.7%), but its shape is still often
preferred for solvating a single roughly-spherical solute.

Use as `PackmolSystem(structure_types; unitcell=octahedral_unitcell(Float64, d), unitcell_center=[cx,cy,cz], ...)`
— equivalent to the input file's `pbc octahedral cx cy cz d`.
"""
function octahedral_unitcell(::Type{T}, d) where {T}
    α = acosd(one(T) / 3)
    β = acosd(-one(T) / 3)
    Matrix{T}(unitcell_matrix(T, d, d, d, α, β, α))
end

#=
    parse_pbc_octahedral(::Type{T}, values, D::Int) where {T}

Parses the `cx cy cz d` tail of a `pbc octahedral cx cy cz d` input file
line into `(unitcell, unitcell_center)`, as consumed by `read_packmol_input`.
=#
function parse_pbc_octahedral(::Type{T}, values, D::Int) where {T}
    D == 3 || throw(ArgumentError("pbc octahedral is only supported for D=3, got D=$D"))
    length(values) == 4 || throw(ArgumentError(
        "pbc octahedral requires 4 values (cx cy cz d), got $(length(values))"
    ))
    cx, cy, cz, d = (_parse_value(T, "pbc", v) for v in values)
    unitcell = octahedral_unitcell(T, d)
    unitcell_center = SVector{D,T}(cx, cy, cz)
    return unitcell, unitcell_center
end

"""
    triclinic_to_octahedral(x, unitcell, center)
    triclinic_to_octahedral(x, packmol_system::PackmolSystem)

Remaps position(s) `x` to the periodic image closest to `center` — the
truncated-octahedron-shaped (Wigner-Seitz) fundamental domain of the cell,
as opposed to the skewed-parallelepiped one `wrap_to_center` produces. `x`
can be a single `SVector{3}` or a vector of them. `unitcell` and `center`
can also be given together as a `PackmolSystem`.

See also [`octahedral_to_triclinic`](@ref) for the reverse mapping.
"""
triclinic_to_octahedral(x::SVector{3,T}, unitcell::AbstractMatrix, center::SVector{3,T}) where {T} =
    _nearest_periodic_image(x, unitcell, center)

triclinic_to_octahedral(xs::AbstractVector{<:SVector{3}}, unitcell::AbstractMatrix, center::SVector{3}) =
    [triclinic_to_octahedral(x, unitcell, center) for x in xs]

triclinic_to_octahedral(x, packmol_system::PackmolSystem) =
    triclinic_to_octahedral(x, packmol_system.unitcell, packmol_system.unitcell_center)

"""
    triclinic_to_octahedral(packmol_system::PackmolSystem)

Sets `packmol_system.periodic_boundary_style = :octahedral` and returns
`packmol_system` — unlike the coordinate-remapping methods above, this one
touches no coordinates itself. It just selects, for subsequent calls to
[`get_atoms`](@ref)/`write_output`, which fundamental-domain shape molecule
positions get wrapped into: the truncated-octahedron (Wigner-Seitz)
"compact" shape instead of the default triclinic one. The actual remapping
happens on demand each time `get_atoms` runs, not here.

Requires `packmol_system.unitcell` to already be set (typically to an
octahedral cell — e.g. via [`octahedral_unitcell`](@ref) or the input file's
`pbc octahedral ...` — though any triclinic cell works).
"""
function triclinic_to_octahedral(packmol_system::PackmolSystem)
    isnothing(packmol_system.unitcell) && throw(ArgumentError(
        "triclinic_to_octahedral(::PackmolSystem) requires packmol_system.unitcell to be set (no PBC is configured)."
    ))
    packmol_system.periodic_boundary_style = :octahedral
    return packmol_system
end

"""
    octahedral_to_triclinic(x, unitcell, center)
    octahedral_to_triclinic(x, packmol_system::PackmolSystem)

The other direction of [`triclinic_to_octahedral`](@ref): remaps
position(s) `x` into the skewed-parallelepiped fundamental domain of the
cell (the ordinary `unitcell` representation, e.g. what CRYST1/PDB output
implies). This is just `wrap_to_center` under a name that pairs with
`triclinic_to_octahedral` for this box. `x` can be a single `SVector{3}`
or a vector of them; `unitcell` and `center` can also be given together as
a `PackmolSystem`.
"""
octahedral_to_triclinic(x::SVector{3,T}, unitcell::AbstractMatrix, center::SVector{3,T}) where {T} =
    wrap_to_center(x, unitcell, center)

octahedral_to_triclinic(xs::AbstractVector{<:SVector{3}}, unitcell::AbstractMatrix, center::SVector{3}) =
    [octahedral_to_triclinic(x, unitcell, center) for x in xs]

octahedral_to_triclinic(x, packmol_system::PackmolSystem) =
    octahedral_to_triclinic(x, packmol_system.unitcell, packmol_system.unitcell_center)

"""
    octahedral_to_triclinic(packmol_system::PackmolSystem)

The other direction of [`triclinic_to_octahedral`](@ref)`(packmol_system)`:
sets `packmol_system.periodic_boundary_style = :triclinic` (the default) and
returns `packmol_system`, so [`get_atoms`](@ref)/`write_output` go back to
wrapping molecule positions into the ordinary skewed-parallelepiped shape.
"""
function octahedral_to_triclinic(packmol_system::PackmolSystem)
    packmol_system.periodic_boundary_style = :triclinic
    return packmol_system
end
