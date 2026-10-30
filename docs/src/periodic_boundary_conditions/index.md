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
