# Defining systems in Julia

This page documents building and packing a system entirely from Julia code, without a
text input file. See [Input files](../input_files/index.md) for the text-file equivalent, and
[Constraints](../constraints/index.md) for the full list of constraint constructors
(`InsideBox`, `OutsideSphere`, `InsideCylinder`, ...) used below.

A system is built from the bottom up: one [`structure_type`](@ref) per kind
of molecule, collected into a [`PackmolSystem`](@ref), then packed with
[`packmol`](@ref).

```julia
using Packmol

water = structure_type(
    "water.pdb";
    number = 1000,
    constraints = [InsideBox([0.,0.,0.],[40.,40.,40.])],
)

sys = PackmolSystem([water]; output = "box.pdb", tolerance = 2.0, seed = 1234)

packmol(sys)
```

This is equivalent to the input file:

```
tolerance 2.0
output box.pdb
seed 1234

structure water.pdb
    number 1000
    inside box 0. 0. 0. 40. 40. 40.
end structure
```
