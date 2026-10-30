# Box

```
inside box  xmin ymin zmin  xmax ymax zmax
outside box xmin ymin zmin  xmax ymax zmax
```

An axis-aligned rectangular box spanning the given corners.

## Example

300 water molecules (about 1 g/mL) in a box with sides of 30, 20 and 15 Å:

```
structure water.pdb
    number 300
    inside box 0. 0. 0. 30. 20. 15.
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl):

```@setup box
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example box
using Packmol, PDBTools
water = structure_type("water.pdb";
    number = 300,
    constraints = [InsideBox([0., 0., 0.], [30., 20., 15.])],
)
system = PackmolSystem([water]; tolerance = 2.0, output = "box.pdb")
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system); unitcell = UnitCell([30., 20., 15.]))
```
