# Global keywords

| Keyword | Meaning |
|:--|:--|
| `tolerance <value>` | Minimum distance between atoms of different molecules (required). |
| `output <file>` | Path of the output structure file. |
| `filetype pdb\|xyz` | Format of the structure and output files (default `pdb`). Files with the extensions `.pdb`, `.xyz`, `.cif` or `.mmcif` (mmCIF, structure files only) are read and written in that format, whatever the `filetype`. |
| `seed <int>` | Random number generator seed (default `1234567`). |
| `movefrac <value>` | Maximum fraction of the molecules of each type moved by the bad-molecule heuristic in each loop (default `0.05`). The `movefrac` keyword argument of [`packmol`](@ref) takes precedence. |
| `nthreads <n>` | Maximum number of threads used (default: all). The `nthreads` keyword argument of [`packmol`](@ref) and the `-t` flag of the executable take precedence. |
| `nloop <int>` | Maximum number of packing loops (default `200`). The `nloop` keyword argument of [`packmol`](@ref) takes precedence. |
| `maxit <int>` | Optimizer iterations per packing loop (default `800`). The `maxit` keyword argument of [`packmol`](@ref) takes precedence. |
| `radscale <value>` | Loose-start radius scale factor (see [`packmol`](@ref)); `discale` is accepted as a legacy alias. |
| `tolerance_precision <value>`, `constraint_precision <value>` | Convergence tolerances for the distance and constraint checks. |
| `stall_tolerance <value>` | Relative-improvement threshold below which a packing chunk is considered stalled and cut short (default `1e-2 * tolerance_precision`). |
| `max_random_init <int>` | Best-of-`N` random placements tried per molecule during initialization. |
| `adjust_constraints_on_init yes\|no` | Whether to run constraint-only pre-optimization before packing. |
| `avoid_overlap yes\|no` | Whether trial placements are rejected if they overlap fixed atoms. |
| `connect yes\|no` | Whether to write the connectivity of the structure files (their `CONECT` records) to the output, for every molecule (default `yes`). `ignore_conect` (no value) is the same as `connect no`. |
| `writeout <int>` | Write intermediate output every `N` loops. |
| `writebad` | Flag (no value): also write output if packing did not converge. |
| `randominitialpoint` | Flag: force a fully random initial point instead of the direct-sampling shortcut. |
| `add_amber_ter` | Flag: add a `TER` record after every molecule (AMBER convention). |
| `amber_ter_preserve` | Flag: preserve the `TER` records of the structure files, in every copy of the molecule. |
| `add_box_sides` | Flag: add a `CRYST1` record with the box dimensions. Added automatically (without needing this flag) whenever `pbc`/`unitcell` is set — with the real triclinic angles in that case, rather than the extrema-derived, always-orthorhombic box this flag alone falls back to. |
| `optimization_print_level <int>` | Verbosity of the optimizer's own output. |
| `chkgrad` | Flag: check analytical vs. numerical gradients (development use). |
| `check` | Write the initial approximation only, then stop (no packing). |
| `pbc <a b c>` or `pbc <xmin ymin zmin xmax ymax zmax>` | Orthorhombic periodic box: either side lengths (the box from the origin to `(a, b, c)`, as in Fortran Packmol) or explicit min/max corners. |
| `pbc dodecahedral <cx cy cz d>` | Rhombic dodecahedron periodic box of size `d` centered at `(cx,cy,cz)` — stored (and written out) as the equivalent triclinic cell. |
| `pbc octahedral <cx cy cz d>` | Truncated octahedron periodic box of size `d` centered at `(cx,cy,cz)` — likewise stored as the equivalent triclinic cell. |
| `unitcell <a b c α β γ>` | Triclinic periodic box, CRYST1-style (centered at the origin). |
| `restart_from <file>` | Skip initial placement for every free molecule and read their positions from `<file>` instead — see [Restarting a run](@ref) below. Also valid inside a `structure ... end structure` block, to restart only that structure type. |
| `restart_to <file>` | Write the coordinates of the system to `<file>` (PDB, or XYZ by the extension) each time output is written (including intermediate `writeout` writes). Also valid per structure type. |

See [Periodic boundary conditions](../periodic_boundary_conditions/index.md) for the
full `pbc`/`unitcell` syntax (including the dodecahedral/octahedral shapes)
and how they interact with the [Constraints](../constraints/index.md) declared
inside `structure` blocks.

Fortran Packmol keywords that have no use in Packmol.jl (`fscale`, `sidemax`,
`maxmove`, `nloop0`, `movebadrandom`, `use_short_tol` and the other short-tolerance keywords,
`hexadecimal_indices`, `non_standard_conect`, ...) are accepted with a
deprecation warning and ignored, so that old input files still run. The other
file formats of Fortran Packmol (`filetype tinker`, `filetype moldy`) are not
supported. See
[`dev/drop-in-replacement.md`](https://github.com/m3g/Packmol.jl/blob/main/dev/drop-in-replacement.md)
for the current status.
