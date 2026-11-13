# constraints

The geometric constraints of molecules or of some of their atoms. Each constraint type
defines its penalty and gradient as functions of the position of an atom, and how it is
read from an input file:

- `constraint_penalty(c, x)` and `constraint_gradient(c, x)`: the penalty is quadratic,
  `weight * violation^2`, with the violation about a distance in Å (gradient of norm ≈ 1).
  Never without units, as Fortran's ellipsoid (`ρ² - 1`): for large regions the
  gradient of small violations falls below the stopping tolerance of SPGBox (1e-5) while
  the fits of the initial placement ask for penalties below 1e-10, and they stall. Nor in
  Å², as Fortran's spheres and the radius of its cylinders (`d² - r²`): the penalty grows
  as d⁴, and atoms far out of a large region give huge objectives and steps;
- `parse_constraint(::Val{Symbol("inside sphere")}, data, T, filename)`: the constraint of
  an input file line.

| File | Contents |
|---|---|
| `constraints_base.jl` | The abstract `Constraint`, the `Inside`/`Outside`/`Over`/`Below` markers, `show`, the default weights of the penalties (`weight_default`), the generic `parse_constraint`, and the fallback of `_inside_box` (the bounding box of the region of an inside constraint, defined in the file of each constraint with a bounded region; used by the initial placement, `free_volume` and `voxel_place`). |
| `boxes.jl`, `spheres.jl`, `planes.jl`, `cylinders.jl`, `ellipsoids.jl` | The constraints of Fortran Packmol. |
| `gaussians.jl`, `waves.jl`, `exponentials.jl` | Curved surfaces (Packmol.jl only), generalizing the planes, and their radial versions. |
| `surfaces.jl` | Surfaces (or curves) given by points, with or without normals (`PointSurface`, `read_surface`, TS2CG point files), the `within surface`/`not within surface` constraints of the distance to them (`SurfaceDistance`), and the `inside surface`/`outside surface` constraints of the sides of a closed surface with outward normals (`SurfaceSide`): a smooth signed distance (implicit moving least squares) near the points, cell lists of the points, a descent on the points for the nearest one far from them, and the side of each empty cell of the grid, so that positions far from the surface cost a lookup. Not isbits (it holds the points): the union is not stored inline. |
| `custom.jl` | `CustomConstraint`: the region `f(x) ≤ 0` of a function of the user (Julia only, not in the executable). |
| `constraint_types.jl` | `AnyConstraint{T}`, the closed union of the constraint types (the element type of the constraint lists), and `_constraint_dispatch`, through which functions are called on its elements (the union is too large for automatic union splitting). |
| `rotations.jl` | Rotation constraints (`constrain_rotation`): the angle between an axis and a vector of the molecule. |

A new constraint type needs: the type and its penalty and gradient, `parse_constraint`
methods for its keywords, an entry in the union of `constraint_types.jl`, tests of its
gradient (`test/constraints/`), and a documentation page with an example (see `CLAUDE.md`).
