# Rhombic dodecahedral box

`pbc dodecahedral cx cy cz d` sets up a rhombic dodecahedron cell of "size"
`d` centered at `(cx,cy,cz)` — internally just the triclinic cell `a = b =
c = d`, `α = β = 60°`, `γ = 90°` (the standard MD convention, matching
GROMACS's `editconf -bt dodecahedron -d d`), so it's written out (e.g. in
`CRYST1`) and behaves exactly like any other `unitcell`. For a given cell
volume, this shape gets a bigger minimum-image distance than a cube would
(`d` is that guaranteed minimum-image distance, in every direction) — about
70.7% of a cube's volume is enough to enclose the same sphere, which is why
it's the usual choice for solvating a single roughly-spherical solute.

In the Julia API, [`dodecahedral_unitcell`](@ref)`(T, d)` builds the same
3×3 matrix for use as the `unitcell` keyword of [`PackmolSystem`](@ref)
(paired with `unitcell_center = [cx,cy,cz]`). The molecules can be written in the
parallelepiped or in the dodecahedral shape of the cell (see [Output shape](#Output-shape)).

## Example

640 water molecules (about 1 g/mL) in a rhombic dodecahedron with `d = 30` Å:

```
pbc dodecahedral 0. 0. 0. 30.

structure water.pdb
    number 640
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl). By default, the
molecules are written in the parallelepiped (triclinic) shape of the cell, whose edges
are drawn:

```@setup dodecahedral
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example dodecahedral
using Packmol, PDBTools
water = structure_type("water.pdb"; number = 640)
system = PackmolSystem([water]; tolerance = 2.0, output = "dodecahedral.pdb",
    unitcell = dodecahedral_unitcell(Float64, 30.))
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system);
    unitcell = UnitCell(system.unitcell),
    unitcell_origin = -sum(eachcol(system.unitcell)) / 2,
)
```

## Output shape

The parallelepiped and the rhombic dodecahedron are two choices of the fundamental domain of the
same periodic cell: same volume and same system, with a different boundary into which
the molecules are wrapped. The rhombic dodecahedron is the Wigner-Seitz cell: each molecule is at
the periodic image closest to the center of the cell, so the system has the more
compact shape of the rhombic dodecahedron, which is often easier to look at and analyze
(analogous to GROMACS's `trjconv -ur compact`).

[`triclinic_to_dodecahedral`](@ref)`(system)` selects that shape for the output of the
system: [`get_atoms`](@ref) and [`write_output`](@ref) then wrap the molecules into it.
It changes no coordinates and does not affect the packing (see [Relationship with
constraints](relationship_with_constraints.md)); only the output is affected:

```@example dodecahedral
triclinic_to_dodecahedral(system)
visualize(get_atoms(system))
```

[`dodecahedral_to_triclinic`](@ref)`(system)` goes back to the parallelepiped shape (the
default). Both functions also convert coordinates directly:
`triclinic_to_dodecahedral(x, unitcell, center)` and `dodecahedral_to_triclinic(x, unitcell, center)`
map a position, or a vector of positions, `x` between the two shapes.

In an input file, the output shape is chosen with the `periodic_boundary_style` keyword:

```
pbc dodecahedral 0. 0. 0. 30.
periodic_boundary_style dodecahedral
```
