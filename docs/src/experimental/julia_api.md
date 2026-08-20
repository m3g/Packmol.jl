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

## Rotation constraints

[`structure_type`](@ref)'s `constrain_rotation` keyword bounds a molecule's
rotation about one or more axes to `center ± halfwidth` degrees — a hard
bound on the optimizer's own rotation-angle variable, not a soft penalty
(matching the input file's `constrain_rotation <axis> <center_deg>
<halfwidth_deg>`, one call per axis instead of one line per axis):

```julia
lipid = structure_type(
    "lipid.pdb";
    number = 100,
    constrain_rotation = Dict(:z => (0.0, 15.0)),  # keep roughly upright
    constraints = [InsideBox([0.,0.,0.],[80.,80.,40.])],
)
```

## Restart files

`restart_from`/`restart_to`, on both `PackmolSystem` (whole system) and
`structure_type` (one structure type only), work as in the input file's
[Restarting a run](input_files.md#Restarting-a-run) section — a file path,
either Packmol's own raw restart format or a `.pdb` (dispatched by
extension), skips or shortcuts the initial placement step:

```julia
sys = PackmolSystem(structure_types; output="box.pdb", tolerance=2.0, restart_to="box.restart")
packmol(sys)  # ... later, resume packing (or extend a run that stopped early):
sys2 = PackmolSystem(structure_types; output="box2.pdb", tolerance=2.0, restart_from="box.restart")
packmol(sys2)
```

The Julia API additionally accepts `restart_from` as an already-loaded
`Vector{<:Atom}` (from `PDBTools.read_pdb`, or built up any other way)
instead of a file path, skipping the file entirely:

```julia
using PDBTools: read_pdb
atoms = read_pdb("some_configuration.pdb")
sys3 = PackmolSystem(structure_types; output="box3.pdb", tolerance=2.0, restart_from=atoms)
```
