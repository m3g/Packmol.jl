# Custom constraints

```julia
CustomConstraint(f; gradient = nothing, weight = 5.0)
```

A region defined by a function of the position: the atoms are kept where `f(x) ≤ 0`, for
`x` the position of an atom (an `SVector{3}`, in Å). Any region that can be written as an
inequality can be used, combined with the other constraints of the molecule as usual. This
constraint is available only in Julia code, not in input files or in the `packmol`
executable.

Outside the region, the penalty is `weight * f(x)^2`, so `f` should grow with the
distance to the region. A signed distance is the best choice (as `f(x) = norm(x) - R`, for
a sphere of radius `R`), but any function that is positive outside the region, and grows
smoothly away from it, works. The packing also needs the gradient of `f`:

- if `gradient` is not given, it is computed by central finite differences (six
  evaluations of `f` for each atom outside the region);
- it can be computed exactly by automatic differentiation, for example with
  [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl):
  `gradient = x -> ForwardDiff.gradient(f, x)`;
- or it can be written by hand, for a function with a simple gradient.

```@docs
CustomConstraint
```

## Example

Water in a torus around the `z` axis, with radii of 14 and 6 Å (the distance between the
center of the tube and the axis, and the radius of the tube). The distance of a point to
the circle of the center of the tube, minus the radius of the tube, is the signed distance
to the surface of the torus. The torus has a volume of about 9950 Å³, which holds 330
water molecules at 1 g/mL:

```@setup custom
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example custom
using Packmol, PDBTools, ForwardDiff
torus(x) = hypot(hypot(x[1], x[2]) - 14, x[3]) - 6
water = structure_type("water.pdb";
    number = 330,
    constraints = [CustomConstraint(torus; gradient = x -> ForwardDiff.gradient(torus, x))],
)
system = PackmolSystem([water]; tolerance = 2.0, output = "torus.pdb")
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system))
```

Custom constraints can also be applied to some atoms of each molecule, with
`atom_constraints`, as the other constraints (see [Constraints on some
atoms](../julia_api/constraints_on_some_atoms.md)).
