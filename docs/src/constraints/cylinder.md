# Cylinder

```
inside cylinder  cx cy cz  vx vy vz  radius length
outside cylinder cx cy cz  vx vy vz  radius length
```

A finite, capped cylinder: `(cx,cy,cz)` is the center of one end cap,
`(vx,vy,vz)` is the axis direction (any nonzero vector — it does not need to
be normalized), `radius` is the cylinder's radius, and `length` is how far
the cylinder extends from the given end, along the axis, to the other cap.

## Example

134 water molecules (about 1 g/mL) in a cylinder along x, with a radius of 8 Å and a length of 20 Å:

```
structure water.pdb
    number 134
    inside cylinder 0. 0. 0. 1. 0. 0. 8. 20.
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl):

```@setup cylinder
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example cylinder
using Packmol, PDBTools
water = structure_type("water.pdb";
    number = 134,
    constraints = [InsideCylinder([0., 0., 0.], [1., 0., 0.], 8., 20.)],
)
system = PackmolSystem([water]; tolerance = 2.0, output = "cylinder.pdb")
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system))
```
