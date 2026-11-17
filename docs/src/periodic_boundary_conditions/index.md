# Periodic boundary conditions

This page documents the periodic boundary conditions, shared between
[input files](../input_files/index.md) and the [Julia API](../julia_api/index.md). They follow the
keywords of Fortran [Packmol](http://github.com/m3g/packmol), and add new cell shapes.

```
pbc a b c
pbc xmin ymin zmin  xmax ymax zmax
pbc dodecahedral cx cy cz d
pbc octahedral cx cy cz d
unitcell a b c alpha beta gamma
```

`pbc`/`unitcell` are not shape constraints declared inside a `structure`
block (those are documented on the [Constraints](../constraints/index.md) page) —
they are global keywords, given once per input file (or once per
`PackmolSystem`, via the `unitcell`/`unitcell_center` keywords in the [Julia
API](../julia_api/index.md)), that set up a periodic simulation cell. Every form
above ends up stored the same way: as a 3×3 `unitcell` matrix (its columns
are the cell vectors) plus a `unitcell_center`. `pbc` is a convenience for
the two most common shapes (an orthorhombic box, or a rhombic dodecahedron/
truncated octahedron computed from a single size parameter); `unitcell`
gives the fully general triclinic cell directly.

## [Cell shapes](@id cell_shapes)

The shapes of the periodic cells, and how each is given in an input file, with
[`PackmolSystem`](@ref) (its `unitcell` keyword, a matrix of the cell vectors as columns,
centered at the origin unless `unitcell_center` is given), and in the solvation
[Recipes](../recipes/index.md) (their `pbc` keyword):

| Shape | Input file | `PackmolSystem` | Recipes |
|:--|:--|:--|:--|
| [Cube](orthorhombic_box.md) | `pbc a a a` | `unitcell = [a 0 0; 0 a 0; 0 0 a]` | `pbc = :cubic` (the default) |
| [Orthorhombic box](orthorhombic_box.md) | `pbc a b c`, or `pbc xmin ymin zmin xmax ymax zmax` | `unitcell = [a 0 0; 0 b 0; 0 0 c]` | `pbc = :orthorhombic` |
| [Triclinic cell](triclinic_box.md) | `unitcell a b c alpha beta gamma` | `unitcell = UnitCell(a, b, c, α, β, γ)` | — |
| [Rhombic dodecahedron](rhombic_dodecahedral_box.md) | `pbc dodecahedral cx cy cz d` | `unitcell = dodecahedral_unitcell(Float64, d)` | `pbc = :dodecahedral` |
| [Truncated octahedron](truncated_octahedral_box.md) | `pbc octahedral cx cy cz d` | `unitcell = octahedral_unitcell(Float64, d)` | `pbc = :octahedral` |

`UnitCell` is from [PDBTools.jl](https://github.com/m3g/PDBTools.jl). [`packmol`](@ref) also
takes the `pbc` keyword of the input files, with its values, replacing the cell of the system:
`packmol(system; pbc = [40, 40, 40])`, `packmol(system; pbc = (:octahedral, 0, 0, 0, 80))`.

In the recipes, the size of the cell follows from the solute and the `margin` (or
`box_sides`): the cube has the largest side of the box of the solute plus the margins, the
orthorhombic box each side computed independently, and the rhombic dodecahedron and the
truncated octahedron are the smallest that keep every atom of the solute at least `margin`
from their faces. These two are closer to spheres than the cube, and solvate a single,
roughly spherical, solute with fewer solvent molecules. `Membrane` and `MembraneProtein` are
always orthorhombic.
