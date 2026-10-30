```@meta
CollapsedDocStrings = true
```
# Running Packmol

A recipe can be packed directly, without ever writing an explicit input
file:

```julia
using Packmol
sys = SolutionBoxUS(...)
packmol(sys; margin=20.0, output="system.pdb")
```

This builds a `PackmolSystem` directly in memory (the same structure the
native, pure-Julia [Julia API](../julia_api/index.md) uses) and hands it to
`packmol(::PackmolSystem)` — no `.inp` file is ever written; `output` is the
only file produced. Keyword arguments specific to each recipe
(`margin`/`box_sides`/`pbc`, plus `concentration`/`ionic_concentration`
where relevant) are documented under each recipe below; any other keyword
(`nloop`, `iprint`, `seed`, ...) is forwarded to the packing
engine itself — see `packmol(::PackmolSystem)`.

`pbc` selects the shape of the periodic cell that `margin`/`box_sides` size:
`:cubic` (the default) forces all 3 sides to their maximum, giving a cube;
`:orthorhombic` keeps the (possibly unequal) sides as computed; and
`:dodecahedral`/`:octahedral` build a rhombic dodecahedron/truncated
octahedron cell of that same maximum side (see
[`dodecahedral_unitcell`](@ref)/[`octahedral_unitcell`](@ref) and [Periodic
boundary conditions](../periodic_boundary_conditions/index.md)) — smaller-volume,
more sphere-like alternatives to the cube, useful for solvating a single
roughly-spherical solute with fewer solvent molecules. In every case the
cell's volume is `det(unitcell)` (rather than simply the product of 3 sides,
once the shape isn't an orthorhombic box), and that volume is what's used to
size the number of solvent/water/ion molecules needed.

Alternatively, `write_packmol_input` generates the `.inp` file on its own,
to be run (or inspected, or edited) separately — either with the native
engine:

```julia
using Packmol
packmol("box.inp")
```

or with the legacy Fortran binary, wrapped by [`run_packmol_legacy`](@ref):

```julia
using Packmol, Packmol_jll
run_packmol_legacy("box.inp")
```
