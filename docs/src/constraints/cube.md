# Cube

```
inside cube  x y z  side
outside cube x y z  side
```

An axis-aligned cube of the given `side` length, with `(x,y,z)` as its
lower corner (matching Fortran Packmol; note this differs from `box`, whose
first three numbers are also a corner but combined with a second corner
rather than a side length).

## Example

200 water molecules in a cube with sides of 20 Å, with a corner at the origin:

```
structure water.pdb
    number 200
    inside cube 0. 0. 0. 20. # adjust for correct water density
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl):

```@setup cube
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example cube
using Packmol, PDBTools
water = structure_type("water.pdb";
    number = 200,
    constraints = [InsideCube([0., 0., 0.], 20.)],
)
system = PackmolSystem([water]; tolerance = 2.0, output = "cube.pdb")
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system); unitcell = UnitCell([20., 20., 20.]))
```
