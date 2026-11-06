# Global keywords

| Keyword | Meaning |
|:--|:--|
| `tolerance <value>` | Minimum distance between atoms of different molecules (required). |
| `output <file>` | Path of the output structure file. |
| `filetype pdb\|xyz\|mmcif` | Format of the structure and output files (default `pdb`). Files with the extensions `.pdb`, `.xyz`, `.cif` or `.mmcif` are read and written in that format, whatever the `filetype`: `output system.cif` writes the output as mmCIF, with the periodic cell in its `_cell` fields. |
| `seed <int>` | Random number generator seed (default `1234567`); `-1` takes a random seed at each run (see [Replicas of a system](@ref)). |
| `movefrac <value>` | Maximum fraction of the molecules of each type moved by the bad-molecule heuristic in each loop (default `0.05`). The `movefrac` keyword argument of [`packmol`](@ref) takes precedence. |
| `periodic_boundary_style <triclinic\|dodecahedral\|octahedral>` | Shape in which the molecules are written, with periodic boundary conditions: the parallelepiped of the cell (`triclinic`, the default), or the compact rhombic dodecahedron or truncated octahedron (see the [Output shape](../periodic_boundary_conditions/rhombic_dodecahedral_box.md#Output-shape)). |
| `nthreads <n>` | Maximum number of threads used (default: all). The `nthreads` keyword argument of [`packmol`](@ref) and the `-t` flag of the executable take precedence. |
| `nloop <int>` | Maximum number of packing loops (default `500`). The `nloop` keyword argument of [`packmol`](@ref) takes precedence. |
| `maxit <int>` | Optimizer iterations per packing loop (default `30`; short loops, in which the overlapping molecules are moved often, pack fastest). The `maxit` keyword argument of [`packmol`](@ref) takes precedence. |
| `radscale <value>` | Loose-start radius scale factor (see [`packmol`](@ref)): the radii of the atoms are multiplied by it at the start of the packing, decaying to their values in `radscale_decay_loops` loops. If not set, each structure type has its own: `5.0` for free molecules of more than 20 atoms (lipids, polymers, liquid crystals, which otherwise interpenetrate and pack slowly), `1.0` for smaller molecules (water, ions, small solvents, which pack fastest at their radii) and for fixed molecules. If set, it applies to all the structure types without a `radscale` of their own (see [Structure blocks](@ref)). `discale` is accepted as a legacy alias. |
| `writeout <seconds>` | Minimum time between two writes of the output file during the packing (default `60`): the best solution so far is written after a loop that improves it, if the last write ended at least this long before. The final solution is always written. The keyword argument of [`packmol`](@ref) takes precedence. |
| `radscale_decay_loops <int>` | Number of packing loops in which the radii decay from `radscale` times their values to their values, exponentially (default: `20`). The keyword argument of [`packmol`](@ref) takes precedence. |
| `tolerance_precision <value>`, `constraint_precision <value>` | Convergence tolerances for the distance and constraint checks. |
| `max_random_init <int>` | Best-of-`N` random placements tried per molecule during initialization. |
| `adjust_constraints_on_init yes\|no` | Whether to run constraint-only pre-optimization before packing. |
| `avoid_overlap yes\|no` | Whether the initial placement avoids overlaps with the fixed molecules (default `yes`). The molecules moved during the packing always avoid them. |
| `connect yes\|no` | Whether to write the connectivity of the structure files (their `CONECT` records) to the output, for every molecule (default `yes`). `ignore_conect` (no value) is the same as `connect no`. |
| `add_amber_ter` | Flag: add a `TER` record after every molecule (AMBER convention). |
| `amber_ter_preserve` | Flag: preserve the `TER` records of the structure files, in every copy of the molecule. |
| `add_box_sides` | Flag: add a `CRYST1` record with the box dimensions. Added automatically (without needing this flag) whenever `pbc`/`unitcell` is set — with the real triclinic angles in that case, rather than the extrema-derived, always-orthorhombic box this flag alone falls back to. |
| `check_initial_placement` | Flag: write the initial placement of the molecules (the initial approximation) to the output file, then stop, without packing: to check that the constraints place the molecules where intended. `check`, Fortran's name, is accepted with a warning. |
| `continue_failed_initial_packing` | Flag: pack even if the initial approximation does not satisfy the constraints. By default, Packmol stops there, writes the initial approximation to the output file and reports the error (exit code 171): constraints that cannot be satisfied even without the distances between the molecules are usually wrong. |
| `pbc <a b c>` or `pbc <xmin ymin zmin xmax ymax zmax>` | Orthorhombic periodic box: either side lengths (the box from the origin to `(a, b, c)`, as in Fortran Packmol) or explicit min/max corners. |
| `pbc dodecahedral <cx cy cz d>` | Rhombic dodecahedron periodic box of size `d` centered at `(cx,cy,cz)` — stored (and written out) as the equivalent triclinic cell. |
| `pbc octahedral <cx cy cz d>` | Truncated octahedron periodic box of size `d` centered at `(cx,cy,cz)` — likewise stored as the equivalent triclinic cell. |
| `unitcell <a b c α β γ>` | Triclinic periodic box, CRYST1-style (centered at the origin). |
| `restart_from <file>` | Skip initial placement for every free molecule and read their positions from `<file>` instead — see [Restarting a run](@ref) below. Also valid inside a `structure ... end structure` block, to restart only that structure type. |

See [Periodic boundary conditions](../periodic_boundary_conditions/index.md) for the
full `pbc`/`unitcell` syntax (including the dodecahedral/octahedral shapes)
and how they interact with the [Constraints](../constraints/index.md) declared
inside `structure` blocks.

Fortran Packmol keywords that have no use in Packmol.jl (`fscale`, `sidemax`,
`maxmove`, `nloop0`, `movebadrandom`, `use_short_tol` and the other short-tolerance keywords,
`hexadecimal_indices`, `non_standard_conect`, `randominitialpoint`,
`optimization_print_level`, `chkgrad`, `restart_to`, ...) are accepted with a
deprecation warning and ignored, so that old input files still run. `writebad` is
ignored too: the best solution is written when the packing does not converge. The other
file formats of Fortran Packmol (`filetype tinker`, `filetype moldy`) are not
supported. See
[`dev/drop-in-replacement.md`](https://github.com/m3g/Packmol.jl/blob/main/dev/drop-in-replacement.md)
for the current status.
