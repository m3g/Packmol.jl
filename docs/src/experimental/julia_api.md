# Defining systems in Julia

!!! warning "Experimental"
    This documents building and packing a system entirely from Julia code,
    without a text input file. See [Input files](input_files.md) for the
    text-file equivalent, and [Constraints](constraints.md) for the full
    list of constraint constructors (`InsideBox`, `OutsideSphere`,
    `InsideCylinder`, ...) used below.

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

## Multiple structure types

Pass one `structure_type` per kind of molecule; each gets its own
constraints:

```julia
solute = structure_type("solute.pdb"; number=1, fixed=([0.,0.,0.],[0.,0.,0.]))
water = structure_type(
    "water.pdb";
    number = 1000,
    constraints = [OutsideSphere([0.,0.,0.],15.), InsideBox([-30.,-30.,-30.],[30.,30.,30.])],
)

sys = PackmolSystem([solute, water]; output = "solvated.pdb", tolerance = 2.0)
packmol(sys)
```

## Combining constraints

Any number of constraints can be given for a `structure_type`; a molecule
must satisfy all of them simultaneously. `InsideBox`/`OutsideSphere`/... are
regular Julia values, so they can be built programmatically (in a loop, from
a function, ...) instead of being written out by hand:

```julia
shells = [
    (r_in, r_out) -> structure_type(
        "water.pdb";
        number = 500,
        constraints = [InsideSphere([0.,0.,0.],r_out), OutsideSphere([0.,0.,0.],r_in)],
    )
    for (r_in, r_out) in [(0.,120.), (150.,200.)]
]
```

## Global options

Any [`PackmolSystem`](@ref) field can be passed as a keyword to the
`PackmolSystem` constructor — `tolerance`, `seed`, `radscale`,
`tolerance_precision`, `unitcell`, and so on, matching the input file's
[Global keywords](@ref) of the same name. [`packmol`](@ref) itself also takes
a few keywords directly (`nloop`, `maxit`, `iprint`, `parallel`,
...) that don't have an input-file equivalent, since they control this
particular run rather than being system properties.
