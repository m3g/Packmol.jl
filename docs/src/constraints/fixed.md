# Fixed molecules

```
fixed x y z  beta gamma theta
center
centerofmass
```

A structure with `fixed` is not packed: it is placed once, at the position `(x, y, z)`
and with the rotation given by the Euler angles `beta`, `gamma` and `theta` (in radians,
as in Fortran Packmol), and the other molecules are packed around it. `number` is
optional (a fixed structure is a single molecule).

The rotation is applied to the coordinates of the structure file, about the origin of
those coordinates, and then the molecule is translated by `(x, y, z)`. With `center` the
molecule is first centered at its geometric center, and with `centerofmass` at its center
of mass (with the masses of the atoms; Fortran Packmol uses unit masses, that is, the
geometric center), so that this center ends up at `(x, y, z)`. `fixed 0. 0. 0. 0. 0. 0.`
without `center` keeps the molecule exactly where it is in its file.

A fixed molecule takes no constraints. With [periodic boundary
conditions](../periodic_boundary_conditions/index.md), its atoms must be inside the
periodic cell.

In the [Julia API](../julia_api/index.md), `fixed = (position, angles)` and
`center = :geometric` or `center = :mass` are keyword arguments of
[`structure_type`](@ref).

## Example

A protein fixed with its geometric center at the origin, surrounded by 1000 water
molecules in a sphere with a radius of 30 Å:

```
structure protein.pdb
    number 1
    center
    fixed 0. 0. 0. 0. 0. 0.
end structure

structure water.pdb
    number 1000
    inside sphere 0. 0. 0. 30.
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl):

```@setup fixed
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "protein.pdb"), "protein.pdb"; force = true)
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example fixed
using Packmol, PDBTools
protein = structure_type("protein.pdb";
    number = 1,
    fixed = ([0., 0., 0.], [0., 0., 0.]),
    center = :geometric,
)
water = structure_type("water.pdb";
    number = 1000,
    constraints = [InsideSphere([0., 0., 0.], 30.)],
)
system = PackmolSystem([protein, water]; tolerance = 2.0, output = "fixed.pdb")
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system); water = true)
```
