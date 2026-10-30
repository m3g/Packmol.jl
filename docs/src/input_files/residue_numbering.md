# Residue numbering

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
