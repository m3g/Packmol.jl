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
| `add_box_sides` | Flag: add a `CRYST1` record with the box dimensions. Added automatically (without needing this flag) whenever `pbc`/`unitcell` is set — with the real triclinic angles in that case, rather than the extrema-derived, always-orthorhombic box this flag alone falls back to. |
| `optimization_print_level <int>` | Verbosity of the optimizer's own output. |
| `chkgrad` | Flag: check analytical vs. numerical gradients (development use). |
| `check` | Write the initial approximation only, then stop (no packing). |
| `pbc <a b c>` or `pbc <xmin ymin zmin xmax ymax zmax>` | Orthorhombic periodic box: either side lengths (centered at the origin) or explicit min/max corners. |
| `unitcell <a b c α β γ>` | Triclinic periodic box, CRYST1-style (centered at the origin). |
| `restart_from <file>` | Skip initial placement for every free molecule and read their positions from `<file>` instead — see [Restarting a run](@ref) below. Also valid inside a `structure ... end structure` block, to restart only that structure type. |
| `restart_to <file>` | Write every free molecule's position to `<file>` each time output is written (including intermediate `writeout` writes). Also valid per structure type. |

A handful of Fortran Packmol keywords are recognized but not yet implemented
(`movefrac`, `movebadrandom`, `use_short_tol` and friends, `nloop`, ...) —
using one prints a warning and it is ignored. See the "Missing keyword
parsing" checklist in the repository's
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
    [chain <c> | changechains]
    [constrain_rotation <x|y|z> <center_deg> <halfwidth_deg>]
    [restart_from <file>]
    [restart_to <file>]
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
| `resnumbers <n>` | Residue numbering scheme for the output — see [Residue numbering](@ref) below. |
| `chain <c>` | Fix every molecule of this type to chain identifier `<c>` in the output, instead of the default auto-assigned one. |
| `changechains` | Alternate between two auto-assigned chain identifiers every other molecule, instead of one chain for the whole type. Not compatible with `chain`. |
| `constrain_rotation <axis> <center_deg> <halfwidth_deg>` | Restrict this molecule's rotation about `<axis>` (`x`, `y`, or `z`) to `<center_deg> ± halfwidth_deg` degrees, as a hard bound on the optimizer's own rotation-angle variable — not a soft penalty. Repeat once per constrained axis. |
| `restart_from <file>` | Restart only this structure type's molecules — see [Restarting a run](@ref). |
| `restart_to <file>` | Write only this structure type's molecules to a restart file — see [Restarting a run](@ref). |
| `atoms <indices> ... end` | Restrict the enclosed lines (`radius`, constraints) to the listed 1-based atom indices instead of the whole molecule. |

Constraint lines (`inside box`, `outside sphere`, `above plane`, ...) are
documented on the [Constraints](constraints.md) page — any number of them can
appear in a structure block, optionally inside an `atoms ... end` block.

## Residue numbering

`resnumbers` controls how the output PDB's residue numbers are assigned, per
structure type:

| Mode | Residue number |
|:--|:--|
| `0` | Sequential per molecule, restarting at 1 for each structure type (or a constant `1` for a `fixed` molecule, since there's only one to number). |
| `1` | Copied verbatim from that atom's own residue number in the template PDB — for a multi-residue template (e.g. a protein), every copy keeps the same internal residue numbers. |
| `2` | Like `1`, but offset by a running counter so consecutive molecules of the same type get consecutive, non-overlapping residue numbers instead of repeating them. |
| `3` | Sequential across the *whole system* (all structure types combined), rather than restarting per type like mode `0`. |

All modes wrap at 9999 (the PDB residue-number field's width), matching
Fortran Packmol. If `resnumbers` isn't given, it's chosen automatically —
mode `0` for a single-residue template, mode `1` for a multi-residue one —
matching Fortran's own default.

## Restarting a run

`restart_from`/`restart_to` skip (or shortcut) the initial placement step and
resume packing from a previously-saved configuration — useful to extend a run
that stopped short of `nloop`, or to continue packing after adding more
molecules to the input file. As a global keyword, they cover every free
molecule in the system (and `restart_from` then skips the whole
initial-placement pipeline entirely, not just the random step — the main
reason to use it on a large system); given inside a `structure ... end
structure` block, they cover only that type's molecules instead, with the
rest of the system still placed normally.

`restart_from <file>` accepts two formats, chosen by `<file>`'s extension:

- Anything else: Packmol's own raw restart format (plain text, one line per
  molecule: center of mass and the three Euler angles, in radians) — what
  `restart_to` always writes.
- A `.pdb` file: any PDB with the same number of atoms and molecules as the
  structure type(s) being restarted — including, but not limited to,
  Packmol's own prior output. Each molecule's (center of mass, orientation)
  is recovered by rigid-body-aligning that structure type's own template
  onto the observed atom positions (the Kabsch algorithm), rather than read
  out directly, so the atoms don't need to already sit at exactly the
  template's bond lengths/angles.

Fixed molecules are never included in a restart, in either direction — their
position is already fully determined by their own `fixed` keyword.
