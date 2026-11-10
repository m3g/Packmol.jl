# constraints

The geometric constraints of molecules or of some of their atoms. Each constraint type
defines its penalty and gradient as functions of the position of an atom, and how it is
read from an input file:

- `constraint_penalty(c, x)` and `constraint_gradient(c, x)`: the penalty is quadratic,
  `weight * violation^2`, with the violation about a distance in Å (gradient of norm ≈ 1,
  or larger: the spheres and the radius of the cylinders use `d² - r²`, in Å², as
  Fortran). Never without units, as Fortran's ellipsoid (`ρ² - 1`): for large regions the
  gradient of small violations falls below the stopping tolerance of SPGBox (1e-5) while
  the fits of the initial placement ask for penalties below 1e-10, and they stall;
- `parse_constraint(::Val{Symbol("inside sphere")}, data, T, filename)`: the constraint of
  an input file line.

| File | Contents |
|---|---|
| `constraints_base.jl` | The abstract `Constraint`, the `Inside`/`Outside`/`Over`/`Below` markers, `show`, and the generic `parse_constraint`. |
| `boxes.jl`, `spheres.jl`, `planes.jl`, `cylinders.jl`, `ellipsoids.jl` | The constraints of Fortran Packmol. |
| `gaussians.jl`, `waves.jl`, `exponentials.jl` | Curved surfaces (Packmol.jl only), generalizing the planes, and their radial versions. |
| `custom.jl` | `CustomConstraint`: the region `f(x) ≤ 0` of a function of the user (Julia only, not in the executable). |
| `constraint_types.jl` | `AnyConstraint{T}`, the closed union of the constraint types (the element type of the constraint lists), and `_constraint_dispatch`, through which functions are called on its elements (the union is too large for automatic union splitting). |
| `rotations.jl` | Rotation constraints (`constrain_rotation`): the angle between an axis and a vector of the molecule. |

A new constraint type needs: the type and its penalty and gradient, `parse_constraint`
methods for its keywords, an entry in the union of `constraint_types.jl`, tests of its
gradient (`test/constraints/`), and a documentation page with an example (see `CLAUDE.md`).
