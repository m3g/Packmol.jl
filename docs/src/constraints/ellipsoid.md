# Ellipsoid

```
inside ellipsoid  cx cy cz  a b c  scale
outside ellipsoid cx cy cz  a b c  scale
```

An axis-aligned ellipsoid centered at `(cx,cy,cz)` with semi-axes `a`, `b`,
`c` along x, y, z respectively, uniformly scaled by the dimensionless
`scale` factor — i.e. the effective semi-axes are `a*scale`, `b*scale`,
`c*scale`. `scale` lets the same base ellipsoid (e.g. one fit to a
reference structure) be grown or shrunk without recomputing `a`, `b`, `c`.

## Example

81 water molecules (about 1 g/mL) in an ellipsoid with semi-axes of 12, 8 and 6 Å:

```
structure water.pdb
    number 81
    inside ellipsoid 0. 0. 0. 12. 8. 6. 1.0
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl):

```@setup ellipsoid
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example ellipsoid
using Packmol, PDBTools
water = structure_type("water.pdb";
    number = 81,
    constraints = [InsideEllipsoid([0., 0., 0.], 12., 8., 6., 1.0)],
)
system = PackmolSystem([water]; tolerance = 2.0, output = "ellipsoid.pdb")
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system))
```
