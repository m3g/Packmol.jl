# Plane

```
above plane x y z d
over plane  x y z d
below plane x y z d
```

A half-space cut by the plane `x*a + y*b + z*c = d`, where `(a,b,c)` is the
plane's normal vector (the first three numbers) and `d` is given as the
fourth. `above`/`over` are synonyms (both accepted by the original Fortran
Packmol); `below` is the complementary half-space.

## Example

100 water molecules in the half of a box below the diagonal plane `x + y = 20`:

```
structure water.pdb
    number 100
    inside box 0. 0. 0. 20. 20. 20.
    below plane 1. 1. 0. 20.
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl):

```@setup plane
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example plane
using Packmol, PDBTools
water = structure_type("water.pdb";
    number = 100,
    constraints = [InsideBox([0., 0., 0.], [20., 20., 20.]), BelowPlane([1., 1., 0.], 20.)],
)
system = PackmolSystem([water]; tolerance = 2.0, output = "plane.pdb")
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system); unitcell = UnitCell([20., 20., 20.]))
```
