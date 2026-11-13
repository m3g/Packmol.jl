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

## The scale of `f`

Write `f` in Å, as a distance, not as a ratio without units. The fits of the initial
placement push every molecule until its penalty is practically zero (below 1e-10), and
the optimizer stops when the gradient of the penalty, `2 weight f(x) ∇f(x)`, falls below
1e-5. With `|∇f| ≈ 1` (a distance), that happens at violations of about 1e-6 Å, where the
penalty is negligible. With a function without units, `|∇f|` may be tiny: for
`f(x) = (norm(x)/R)^2 - 1`, the inside of a sphere of radius `R`, it is `2/R`, and with
`R = 320` Å the optimizer stops with the atoms still about 0.01 Å outside, while their
penalties are above 1e-10. The fits then run all their loops without converging, as
packing 20,000 lipids in a shell inside such a sphere shows:

| `f(x)` | units of `f` | `abs(∇f)` at the surface | first fit |
|:-------|:-----:|:-----:|:-----|
| `norm(x) - R` | Å | 1 | converged in 0.4 s |
| `sum(abs2, x) - R^2` | Å² | `2R` | converged in 0.6 s |
| `sum(abs2, x) / R^2 - 1` | none | `2/R` | not converged after 200 loops |

So:

- prefer a signed distance, or a function with a gradient of about 1 near the boundary of
  the region (a function that grows faster, as `sum(abs2, x) - R^2`, also works near the
  region, but far from it the penalty grows as the fourth power of the distance: it
  dominates the distances between atoms in the objective function, and atoms placed
  hundreds of Å out of a large region make the first steps of the optimizer huge);
- do not divide `f` by the size of the region: scale a function given without units by a
  length of the region (the inside of an ellipsoid of semi-axes `a, b, c` is `m (ρ - 1)`,
  for `ρ = norm(x ./ (a, b, c))` and `m = min(a, b, c)`, as the `Ellipsoid` constraint);
- or, which has the same effect, increase `weight` (by the square of that length).

The gradient can be checked with ForwardDiff.jl, comparing
`ForwardDiff.gradient(x -> Packmol.constraint_penalty(c, x), x)` to
`Packmol.constraint_gradient(c, x)` at points outside the region.

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
visualize(get_atoms(system);
    rotate = [(:z, 30), (:x, -70)],
)
```

Custom constraints can also be applied to some atoms of each molecule, with
`atom_constraints`, as the other constraints (see [Constraints on some
atoms](../julia_api/constraints_on_some_atoms.md)).
