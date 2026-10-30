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

## Constraints on some atoms

The `atom_constraints` keyword of [`structure_type`](@ref) applies constraints to some
atoms of the molecule only, as the `atoms ... end atoms` blocks of an input file. It is a
list of `atom indices => constraints` pairs (one constraint or a vector of them), with
the indices of the atoms in the structure file:

```julia
lipid = structure_type("palmitoil.pdb";
    number = 50,
    constraints = [InsideBox([0.,0.,0.],[40.,40.,14.])], # all atoms
    atom_constraints = [
        [31, 32] => BelowPlane([0.,0.,1.], 2.0), # polar head down
        [1, 2] => AbovePlane([0.,0.,1.], 12.0),  # tail up
    ],
)
```

is the input file block:

```
structure palmitoil.pdb
  number 50
  inside box 0. 0. 0. 40. 40. 14.
  atoms 31 32
    below plane 0. 0. 1. 2.
  end atoms
  atoms 1 2
    above plane 0. 0. 1. 12.
  end atoms
end structure
```

See the [Lipid bilayer](../tutorials/bilayer.md) tutorial.

## Global options

Any [`PackmolSystem`](@ref) field can be passed as a keyword to the
`PackmolSystem` constructor — `tolerance`, `seed`, `radscale`,
`tolerance_precision`, `unitcell`, and so on, matching the input file's
[Global keywords](@ref) of the same name. [`packmol`](@ref) itself also takes
a few keywords directly (`nloop`, `maxit`, `iprint`, `parallel`,
...) that don't have an input-file equivalent, since they control this
particular run rather than being system properties.

## Rotation constraints

[`structure_type`](@ref)'s `constrain_rotation` keyword keeps the angle between
an axis and a vector of the molecule within `center ± halfwidth` degrees (see
[Rotation constraint](@ref)). Each axis takes `(center, halfwidth)`, for the
longest axis of the molecule, or `(center, halfwidth, from_atoms, to_atoms)`,
for the vector between the centers of two sets of atoms:

```julia
lipid = structure_type(
    "lipid.pdb";
    number = 100,
    # head (atoms 1 and 2) to tail (atom 42) within 15° of z
    constrain_rotation = Dict(:z => (0.0, 15.0, [1, 2], [42])),
    constraints = [InsideBox([0.,0.,0.],[80.,80.,40.])],
)
```

## Restart files

`restart_from`/`restart_to`, on both `PackmolSystem` (whole system) and
`structure_type` (one structure type only), work as in the input file's
[Restarting a run](input_files.md#Restarting-a-run) section — a coordinate
file (PDB, XYZ or mmCIF, by extension) skips or shortcuts the initial placement
step:

```julia
sys = PackmolSystem(structure_types; output="box.pdb", tolerance=2.0, restart_to="box_restart.pdb")
packmol(sys)  # ... later, resume packing (or extend a run that stopped early):
sys2 = PackmolSystem(structure_types; output="box2.pdb", tolerance=2.0, restart_from="box_restart.pdb")
packmol(sys2)
```

The Julia API additionally accepts `restart_from` as an already-loaded
`Vector{<:Atom}` (from `PDBTools.read_pdb` (PDBTools ≥ 4), or built up any other way)
instead of a file path, skipping the file entirely:

```julia
using PDBTools: read_pdb
atoms = read_pdb("some_configuration.pdb")
sys3 = PackmolSystem(structure_types; output="box3.pdb", tolerance=2.0, restart_from=atoms)
```
