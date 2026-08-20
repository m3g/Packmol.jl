# Input files

!!! warning "Experimental"
    This documents the native Julia packing engine's text input file, read
    with `read_packmol_input` and run with [`packmol`](@ref). It tracks the
    original Fortran [Packmol](http://github.com/m3g/packmol) syntax where a
    keyword is supported, but many Fortran keywords aren't implemented yet
    (noted below). To build a system directly from Julia code instead of a
    text file, see [Defining systems in Julia](julia_api.md).

An input file has two parts: a handful of global keyword lines, and one or
more `structure ... end structure` blocks, one per kind of molecule to pack.
Lines starting with `#` are comments.

```
tolerance 2.0
output box.pdb
seed 1234

structure water.pdb
    number 1000
    inside box 0. 0. 0. 40. 40. 40.
end structure
```

## Global keywords

| Keyword | Meaning |
|:--|:--|
| `tolerance <value>` | Minimum distance between atoms of different molecules (required). |
| `output <file>` | Path of the output structure file. |
| `filetype <type>` | Output file type (default `pdb`). |
| `seed <int>` | Random number generator seed (default `1234567`). |
| `maxit <int>` | Iterations per packing loop chunk. |
| `radscale <value>` | Loose-start radius scale factor (see [`packmol`](@ref)); `discale` is accepted as a legacy alias. |
| `tolerance_precision <value>`, `constraint_precision <value>` | Convergence tolerances for the distance and constraint checks. |
| `stall_tolerance <value>` | Relative-improvement threshold below which a packing chunk is considered stalled and cut short (default `1e-2 * tolerance_precision`). |
| `max_random_init <int>` | Best-of-`N` random placements tried per molecule during initialization. |
| `adjust_constraints_on_init yes\|no` | Whether to run constraint-only pre-optimization before packing. |
| `avoid_overlap yes\|no` | Whether trial placements are rejected if they overlap fixed atoms. |
| `connect yes\|no` | Whether to look for `CONECT` records (not yet implemented — see below). |
| `writeout <int>` | Write intermediate output every `N` loops. |
| `writebad` | Flag (no value): also write output if packing did not converge. |
| `randominitialpoint` | Flag: force a fully random initial point instead of the direct-sampling shortcut. |
| `add_amber_ter yes\|no` | Whether to add AMBER `TER` records. |
| `amber_ter_preserve` | Flag: preserve existing `TER` records. |
| `add_box_sides` | Flag: add a `CRYST1` record with the box dimensions. |
| `optimization_print_level <int>` | Verbosity of the optimizer's own output. |
| `chkgrad` | Flag: check analytical vs. numerical gradients (development use). |
| `check` | Write the initial approximation only, then stop (no packing). |
| `pbc <a b c>` or `pbc <xmin ymin zmin xmax ymax zmax>` | Orthorhombic periodic box: either side lengths (centered at the origin) or explicit min/max corners. |
| `unitcell <a b c α β γ>` | Triclinic periodic box, CRYST1-style (centered at the origin). |

A handful of Fortran Packmol keywords are recognized but not yet implemented
(`movefrac`, `movebadrandom`, `use_short_tol` and friends, `nloop`,
`restart_from`/`restart_to`, ...) — using one prints a warning and it is
ignored. See the "Missing keyword parsing" checklist in the repository's
[`CLAUDE.md`](https://github.com/m3g/Packmol.jl/blob/main/CLAUDE.md) for the
current list.

## Structure blocks

```
structure <pdb file>
    number <n>
    [fixed x y z beta gamma theta]
    [center | centerofmass]
    [radius <value>]
    [resnumbers <0|1|2|3>]
    <constraint lines...>
    [atoms <indices...>
        <per-atom overrides...>
    end]
end structure
```

| Keyword | Meaning |
|:--|:--|
| `number <n>` | Number of copies of this molecule to pack (required). |
| `fixed x y z beta gamma theta` | Don't pack this molecule — place it once at this position and (Euler angle) orientation and leave it fixed. |
| `center` / `centerofmass` | Only valid together with `fixed`: re-center the molecule (geometric center or mass-weighted center of mass) before applying the fixed position/orientation. |
| `radius <value>` | Atom radius for every atom of this structure, overriding the default of `tolerance/2`. |
| `resnumbers <n>` | Residue numbering scheme for the output (parsed, but not yet applied to the output). |
| `atoms <indices> ... end` | Restrict the enclosed lines (`radius`, constraints) to the listed 1-based atom indices instead of the whole molecule. |

Constraint lines (`inside box`, `outside sphere`, `above plane`, ...) are
documented on the [Constraints](constraints.md) page — any number of them can
appear in a structure block, optionally inside an `atoms ... end` block.
