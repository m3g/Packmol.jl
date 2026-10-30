# Choosing the output shape

Called on a `PackmolSystem` directly (no coordinates involved),
`triclinic_to_dodecahedral(packmol_system)`/`triclinic_to_octahedral(packmol_system)`
(and their `..._to_triclinic` inverses) set
`packmol_system.periodic_boundary_style` (`:dodecahedral`/`:octahedral`, or
`:triclinic`, the default) and return `packmol_system` — the shape that
[`get_atoms`](@ref) and [`write_output`](@ref) wrap molecule positions into
when they materialize atomic coordinates, generated fresh from
`packmol_system.molecule_positions` each time rather than stored. No
coordinates are touched by the toggle itself; only the next
`get_atoms`/`write_output` call is affected, and only that call's *output*
— see [Relationship with constraints](@ref) below for why packing itself is
unaffected by this choice.
