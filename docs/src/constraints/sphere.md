# Sphere

```
inside sphere  x y z  radius
outside sphere x y z  radius
```

A sphere of the given `radius` centered at `(x,y,z)`.

## Example

200 water molecules in a sphere with a radius of 12 Å:

```
structure water.pdb
    number 200
    inside sphere 0. 0. 0. 12. # adjust for correct water density
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl):

```@setup sphere
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example sphere
using Packmol, PDBTools
water = structure_type("water.pdb";
    number = 200,
    constraints = [InsideSphere([0., 0., 0.], 12.)],
)
system = PackmolSystem([water]; tolerance = 2.0, output = "sphere.pdb")
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system))
```
