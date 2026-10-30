# Truncated octahedral box

`pbc octahedral cx cy cz d` sets up a truncated octahedron cell of "size"
`d` centered at `(cx,cy,cz)` — again internally just a triclinic cell, `a =
b = c = d`, `α = γ = acosd(1/3) ≈ 70.53°`, `β = acosd(-1/3) ≈ 109.47°`
(matching GROMACS's `editconf -bt octahedron -d d`). It's the Wigner-Seitz
cell of a BCC lattice, as opposed to the dodecahedron's FCC one, and for a
given `d` needs a bit more volume than the dodecahedron does (≈77.0% of a
cube's volume, vs. the dodecahedron's ≈70.7%) — but its shape is still
often preferred for solvating a single roughly-spherical solute.

In the Julia API, [`octahedral_unitcell`](@ref)`(T, d)` builds the same 3×3 matrix
for use as the `unitcell` keyword of [`PackmolSystem`](@ref) (paired with
`unitcell_center = [cx,cy,cz]`). The molecules can be written in the parallelepiped or
in the octahedral shape of the cell (see [Output shape](#Output-shape)).

## Example

690 water molecules (about 1 g/mL) in a truncated octahedron with `d = 30` Å:

```
pbc octahedral 0. 0. 0. 30.

structure water.pdb
    number 690
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl). By default, the
molecules are written in the parallelepiped (triclinic) shape of the cell, whose edges
are drawn:

```@setup octahedral
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example octahedral
using Packmol, PDBTools
water = structure_type("water.pdb"; number = 690)
system = PackmolSystem([water]; tolerance = 2.0, output = "octahedral.pdb",
    unitcell = octahedral_unitcell(Float64, 30.))
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system);
    unitcell = UnitCell(system.unitcell),
    unitcell_origin = -sum(eachcol(system.unitcell)) / 2,
)
```

## Output shape

The parallelepiped and the truncated octahedron are two choices of the fundamental domain of the
same periodic cell: same volume and same system, with a different boundary into which
the molecules are wrapped. The truncated octahedron is the Wigner-Seitz cell: each molecule is at
the periodic image closest to the center of the cell, so the system has the more
compact shape of the truncated octahedron, which is often easier to look at and analyze
(analogous to GROMACS's `trjconv -ur compact`).

[`triclinic_to_octahedral`](@ref)`(system)` selects that shape for the output of the
system: [`get_atoms`](@ref) and [`write_output`](@ref) then wrap the molecules into it.
It changes no coordinates and does not affect the packing (see [Relationship with
constraints](relationship_with_constraints.md)); only the output is affected:

```@example octahedral
triclinic_to_octahedral(system)
visualize(get_atoms(system))
```

[`octahedral_to_triclinic`](@ref)`(system)` goes back to the parallelepiped shape (the
default). Both functions also convert coordinates directly:
`triclinic_to_octahedral(x, unitcell, center)` and `octahedral_to_triclinic(x, unitcell, center)`
map a position, or a vector of positions, `x` between the two shapes.

The output shape is chosen in the Julia API only: the output of an input file is always
in the parallelepiped shape.
