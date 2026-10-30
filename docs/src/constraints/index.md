# Constraints

This page documents the constraint syntax, shared between [input files](../input_files/index.md)
and the [Julia API](../julia_api/index.md). It follows the keywords of Fortran
[Packmol](http://github.com/m3g/packmol), and adds new constraints.

A constraint restricts where the atoms of a structure are allowed to be. In
an input file, each one is declared as a line inside a
`structure ... end structure` block, and takes the form:

```
<placement> <shape> <parameters...>
```

Multiple constraints can be given for the same structure — a molecule must
then satisfy all of them simultaneously (e.g. `inside box ...` together with
`outside sphere ...` carves a box with a spherical hole in it). By default a
constraint applies to every atom of the structure; in an input file, wrapping
a subset of lines in an `atoms <indices> ... end` block restricts a
constraint to just those atoms (see [Structure blocks](../input_files/index.md) in
Input files).

!!! note
    Constraints are independent of, but interact with, periodic boundary
    conditions: setting `pbc`/`unitcell` changes how every constraint below
    is evaluated (atoms are wrapped to the periodic image nearest the cell
    center first), and adds one implicit confinement constraint of its own
    on top of whatever's declared explicitly. See [Periodic boundary
    conditions](../periodic_boundary_conditions/index.md) for the `pbc`/`unitcell`
    keywords themselves and the full relationship between the two.
