# Orthorhombic box

```
pbc a b c
pbc xmin ymin zmin  xmax ymax zmax
```

The `a b c` form is the box of the given side lengths with a corner at the
origin, from `(0, 0, 0)` to `(a, b, c)`, as in Fortran Packmol (it is the same
as `pbc 0 0 0 a b c`). The six-value form gives explicit min/max corners
instead. In both cases the center of the cell is the midpoint of the box,
`((xmin+xmax)/2, ...)`. Either way the resulting cell is diagonal (orthorhombic): no
shear between axes.

## Example

900 water molecules (about 1 g/mL) in a periodic cubic box with sides of 30 Å:

```
pbc 30. 30. 30.

structure water.pdb
    number 900
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl):

```@setup orthorhombic
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example orthorhombic
using Packmol, PDBTools
water = structure_type("water.pdb"; number = 900)
system = PackmolSystem([water]; tolerance = 2.0, output = "orthorhombic.pdb",
    unitcell = [30. 0. 0.; 0. 30. 0.; 0. 0. 30.],
    unitcell_center = [15., 15., 15.])
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system); unitcell = UnitCell(system.unitcell))
```
