# Input files

This page documents the input file of Packmol, read with `read_packmol_input` and run
with [`packmol`](@ref) or the `packmol` executable. It follows the syntax of Fortran
[Packmol](http://github.com/m3g/packmol); the differences are noted below. To build a
system directly from Julia code instead of a text file, see
[Defining systems in Julia](../julia_api/index.md).

An input file has two parts: a handful of global keyword lines, and one or
more `structure ... end structure` blocks, one per kind of molecule to pack.
Anything after a `#` on a line is a comment.

```
tolerance 2.0
output box.pdb
seed 1234

structure water.pdb
    number 1000
    inside box 0. 0. 0. 40. 40. 40.
end structure
```
