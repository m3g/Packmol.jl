```@meta
CollapsedDocStrings = true
```

# Recipes

!!! warning "Experimental"
    Recipes are under active development. Struct fields, keyword arguments,
    and the exact contents of generated input files may still change.

## Overview

A recipe is a higher-level, parameter-driven system setup: instead of
specifying molecule counts directly (as a `structure ... end structure`
block does), a recipe takes target densities or concentrations and works out
the molecule counts and box size itself. Recipes currently cover common
solvation setups:

- [Solute-Solvent system](@ref): a solute solvated by a single solvent.
- [Solute-Solvent-Cossolvent system](@ref): a solute solvated by a mixture of
  two solvents (e.g. water/ethanol), at a target cossolvent concentration.
- [Solute-Water-Ions system](@ref): a solute solvated by water and a
  background salt (e.g. NaCl) at a target ionic concentration, with
  automatic charge neutralization of the solute.
- [Membrane system](@ref): a lipid bilayer or monolayer, optionally mixing
  more than one lipid and more than one solvent, solvated above (and, for a
  bilayer, also below) the membrane.

Vesicles and nanotubes are planned as future recipes; a rhombic
dodecahedron or truncated octahedron box is already available today, for
the three solvation recipes above, via the `pbc` keyword (see [Running
Packmol](@ref)) — `Membrane` is always orthorhombic.

Each solvation recipe is built the same way: construct a `SolutionBoxU*`
data structure describing the components of the system (PDB files, molar
masses, densities), then either call `write_packmol_input` on it to generate
a `.inp` file, or call `packmol` on it directly to build and pack the system
without ever writing one (see [Running Packmol](@ref)).

Concentrations and densities can be given as plain numbers (assumed to be in
the units noted in each function's docstring) or as [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/)
quantities (`55.5u"mol/L"`, `1.0u"g/mL"`, ...). See
[Concentration Unit Conversion](@ref) for the underlying `cconvert` machinery
used to convert between them.

### How to use it

```julia-repl
julia> using Packmol
```

`SolutionBoxUS`, `SolutionBoxUSC`, `SolutionBoxUWI`, and `write_packmol_input`
are exported directly from `Packmol`.

## Running Packmol

A recipe can be packed directly, without ever writing an explicit input
file:

```julia
using Packmol
sys = SolutionBoxUS(...)
packmol(sys; margin=20.0, output="system.pdb")
```

This builds a `PackmolSystem` directly in memory (the same structure the
native, pure-Julia [Julia API](julia_api.md) uses) and hands it to
`packmol(::PackmolSystem)` — no `.inp` file is ever written; `output` is the
only file produced. Keyword arguments specific to each recipe
(`margin`/`box_sides`/`pbc`, plus `concentration`/`ionic_concentration`
where relevant) are documented under each recipe below; any other keyword
(`nloop`, `iprint`, `seed`, ...) is forwarded to the packing
engine itself — see `packmol(::PackmolSystem)`.

`pbc` selects the shape of the periodic cell that `margin`/`box_sides` size:
`:cubic` (the default) forces all 3 sides to their maximum, giving a cube;
`:orthorhombic` keeps the (possibly unequal) sides as computed; and
`:dodecahedral`/`:octahedral` build a rhombic dodecahedron/truncated
octahedron cell of that same maximum side (see
[`dodecahedral_unitcell`](@ref)/[`octahedral_unitcell`](@ref) and [Periodic
boundary conditions](periodic_boundary_conditions.md)) — smaller-volume,
more sphere-like alternatives to the cube, useful for solvating a single
roughly-spherical solute with fewer solvent molecules. In every case the
cell's volume is `det(unitcell)` (rather than simply the product of 3 sides,
once the shape isn't an orthorhombic box), and that volume is what's used to
size the number of solvent/water/ion molecules needed.

Alternatively, `write_packmol_input` generates the `.inp` file on its own,
to be run (or inspected, or edited) separately — either with the native
engine:

```julia
using Packmol
packmol("box.inp")
```

or with the legacy Fortran binary, wrapped by [`run_packmol`](@ref):

```julia
using Packmol
run_packmol("box.inp")
```

## Solute-Solvent system

Here, `SolutionBoxUS` stands for `Solute (U)` and `Solvent (S)`. The density
of the pure solvent is given directly (`g/mL` by default, or `mol/L` for
molarity).

```@docs
SolutionBoxUS
write_packmol_input(::SolutionBoxUS)
packmol(::SolutionBoxUS)
```

### Setting up the system properties

We initialize the system data structure given the PDB files of *one
molecule* of the solute (a polymer, `poly_h.pdb`) and *one molecule* of the
solvent (water):

```@example us
using Packmol
test_dir = Packmol.RecipesDirectory * "/test"
system = SolutionBoxUS(
    solute_pdbfile = "$test_dir/data/poly_h.pdb",
    solvent_pdbfile = "$test_dir/data/water.pdb",
    density = 1.0u"g/mL",
)
```

The molar masses of the solute and solvent can be provided explicitly with
the `solute_molar_mass`/`solvent_molar_mass` keywords. If not, they are
computed from the atom types in the PDB files, which may fail if the mass of
some atom type is unknown.

Finally, we generate an input file for Packmol with:

```@example us
write_packmol_input(
    system;
    margin = 20.0,
    input = "box.inp",
    output = "system.pdb",
)
```

or build and pack the system directly, with `packmol`:

```julia
packmol(system; margin = 20.0, output = "system.pdb")
```

The `input` parameter (`write_packmol_input` only) is the name of the
generated Packmol input file, and `output` is the name assigned to the
packed system.

`margin` sets the size of the box from the solute's own bounding box plus
this margin, in every dimension. By default (`pbc = :cubic`) the box is a
cube (all three sides set to the largest of the three margined dimensions);
`pbc = :orthorhombic` instead sizes each side independently, and
`pbc = :dodecahedral`/`pbc = :octahedral` build a rhombic dodecahedron/
truncated octahedron cell of that same largest side (see [Running
Packmol](@ref) above).

Alternatively, the box size can be given explicitly with
`box_sides = [a, b, c]` (in Å), instead of `margin`.

## Solute-Solvent-Cossolvent system

Here, `SolutionBoxUSC` stands for `Solute (U)`, `Solvent (S)`, and
`Cossolvent (C)`. Concentrations can be given in molarity (`"mol/L"`), molar
fraction (`"x"`), mass fraction (`"w/w"`), or volume fraction (`"v/v"`); see
[Concentration Unit Conversion](@ref) for the full list of unit aliases.

```@docs
SolutionBoxUSC
write_packmol_input(::SolutionBoxUSC)
packmol(::SolutionBoxUSC)
```

### Setting up the system properties

The density of the mixture as a function of cossolvent concentration is
given as a table. Here, for a water/ethanol mixture, as a function of the
molar fraction of ethanol:

```@example usc
using Packmol
density_table = [
#   x cossolvent (ethanol)     density (g/mL)
             0.0000                 0.9981
             0.0416                 0.9820
             0.0890                 0.9685
             0.1434                 0.9537
             0.2066                 0.9369
             0.2809                 0.9151
             0.3695                 0.8923
             0.4769                 0.8685
             0.6098                 0.8450
             0.7786                 0.8195
             1.0000                 0.7906
]
```

We then initialize the system data structure, given the PDB files of *one
molecule* of the solute, the solvent (water), and the cossolvent (ethanol):

```@example usc
test_dir = Packmol.RecipesDirectory * "/test"
system = SolutionBoxUSC(
    solute_pdbfile = "$test_dir/data/poly_h.pdb",
    solvent_pdbfile = "$test_dir/data/water.pdb",
    cossolvent_pdbfile = "$test_dir/data/ethanol.pdb",
    density_table = density_table,
    concentration_units = "x", # molar fraction
)
```

As with `SolutionBoxUS`, molar masses are computed from the PDB files unless
given explicitly.

Finally, we generate the input file at a target cossolvent concentration:

```@example usc
write_packmol_input(
    system;
    concentration = 0.5, # molar fraction of ethanol, by the density_table's default units
    margin = 20.0,
    input = "box.inp",
    output = "system.pdb",
)
```

or, equivalently, pack it directly with `packmol(system; concentration=0.5, margin=20.0, output="system.pdb")`.

`concentration_units` can be passed to interpret `concentration` in units
other than the density table's own (e.g. request a concentration in
`"mol/L"` even though the table above is indexed by molar fraction).
`margin`, `box_sides`, and `pbc` behave as in `SolutionBoxUS`.

## Solute-Water-Ions system

Here, `SolutionBoxUWI` stands for `Solute (U)`, `Water (W)`, and `Ions (I)`:
a solute solvated by water and a background salt at a target ionic
concentration (`0.15u"mol/L"`, physiological saline, by default). The
solute's own net charge (`solute_charge`, computed from the PDB file's
`charge` records unless given explicitly) is automatically neutralized with
extra counter-ions, on top of the bulk salt.

```@docs
SolutionBoxUWI
write_packmol_input(::SolutionBoxUWI)
packmol(::SolutionBoxUWI)
```

### Setting up the system properties

By default the ions are sodium and chloride (charges `+1`/`-1`), generated
automatically, and the solution density as a function of ionic strength is
looked up from a table of aqueous NaCl densities:

```@example uwi
using Packmol
test_dir = Packmol.RecipesDirectory * "/test"
system = SolutionBoxUWI(
    solute_pdbfile = "$test_dir/data/poly_h.pdb",
    solute_charge = 2,
)
```

Custom cation/anion PDB files can be given with the `cation_pdbfile`/
`anion_pdbfile` keywords (with their `cation_charge`/`anion_charge`); a
matching `density_table` should then be provided too, since the default one
is specific to NaCl.

Finally, we generate the input file at a target ionic concentration:

```@example uwi
write_packmol_input(
    system;
    ionic_concentration = 0.15u"mol/L",
    margin = 20.0,
    input = "box.inp",
    output = "system.pdb",
)
```

or, equivalently, pack it directly with `packmol(system; ionic_concentration=0.15u"mol/L", margin=20.0, output="system.pdb")`.

Enough cations or anions (of a single sign) are added on top of the bulk
salt to exactly neutralize `system.solute_charge`; a `SolutionBoxUWI` with
default (monovalent) ions can neutralize any integer solute charge. An
`ArgumentError` is raised if the requested charge cannot be exactly
neutralized with the ions' charges (for instance, a solute charge of `-3`
with only doubly-charged cations available). `margin`, `box_sides`, and
`pbc` behave as in `SolutionBoxUS`.

## Membrane system

`Membrane` builds a lipid bilayer or monolayer under orthorhombic periodic
boundary conditions, with the membrane normal along z. Unlike the solvation
recipes above, every parameter is given up front in the `Membrane`
constructor itself — `write_packmol_input`/`packmol` take only `input`/
`output` (and, for `packmol`, any packing-engine keyword).

```@docs
Membrane
write_packmol_input(::Membrane)
packmol(::Membrane)
```

### Setting up the system properties

Each lipid is identified by a single "head" atom and a single "tail" atom
(1-based indices into its own PDB file); the Euclidean distance between them
sets that lipid's head-to-tail length, and the *longest* one (across all
lipids given) sets the membrane's nominal thickness. Here, a single lipid
type, solvated by water:

```@example membrane
using Packmol
test_dir = Packmol.RecipesDirectory * "/test"
system = Membrane(
    type = :bilayer,
    lipids = "$test_dir/data/lipid.pdb",
    lipid_head = [31],  # the head-group ring
    lipid_tail = [1],   # the far end of the acyl chain
    lipid_molar_ratio = [1.0],
    area_per_lipid = 60.0, # Å², xor total_area
    total_lipids = 200,    # xor total_area, counted over both leaflets
    solvent = "$test_dir/data/water.pdb",
    solvent_layer_width = 20.0, # Å
    solvent_density = 1.0u"g/mL",
)
```

Mixing lipids or solvents just means passing `Vector`s of the same length to
`lipids`/`lipid_head`/`lipid_tail`/`lipid_molar_ratio` (and, respectively,
`solvent`/`solvent_molar_ratio`) instead of scalars/single files — see the
`Membrane` docstring above for the full parameter list, including
`flexibility` (how much tilt/wobble each lipid is allowed while its head and
tail are pinned near their target planes) and the `total_area`/
`total_lipids` either-or.

Finally, generate the input file, or pack it directly:

```@example membrane
write_packmol_input(system; input = "membrane.inp", output = "membrane.pdb")
```

```julia
packmol(system; output = "membrane.pdb")
```
