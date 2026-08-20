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

Membranes, vesicles, nanotubes, and special (octahedral/dodecahedral) box
shapes are planned as future recipes.

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
packmol(sys; margin=20.0, cubic=true, output="system.pdb")
```

This builds a `PackmolSystem` directly in memory (the same structure the
native, pure-Julia [Julia API](julia_api.md) uses) and hands it to
`packmol(::PackmolSystem)` — no `.inp` file is ever written; `output` is the
only file produced. Keyword arguments specific to each recipe
(`margin`/`box_sides`/`cubic`, plus `concentration`/`ionic_concentration`
where relevant) are documented under each recipe below; any other keyword
(`nloop`, `iprint`, `seed`, ...) is forwarded to the packing
engine itself — see `packmol(::PackmolSystem)`.

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
    cubic = true,
    input = "box.inp",
    output = "system.pdb",
)
```

or build and pack the system directly, with `packmol`:

```julia
packmol(system; margin = 20.0, cubic = true, output = "system.pdb")
```

The `input` parameter (`write_packmol_input` only) is the name of the
generated Packmol input file, and `output` is the name assigned to the
packed system.

`margin` sets the size of the box from the solute's own bounding box plus
this margin, in every dimension. If `cubic` is `true` the box is a cube
(all three sides set to the largest of the three margined dimensions);
otherwise it is orthorhombic, with each side sized independently.

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
    cubic = true,
    input = "box.inp",
    output = "system.pdb",
)
```

or, equivalently, pack it directly with `packmol(system; concentration=0.5, margin=20.0, cubic=true, output="system.pdb")`.

`concentration_units` can be passed to interpret `concentration` in units
other than the density table's own (e.g. request a concentration in
`"mol/L"` even though the table above is indexed by molar fraction).
`margin`, `box_sides`, and `cubic` behave as in `SolutionBoxUS`.

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
    cubic = true,
    input = "box.inp",
    output = "system.pdb",
)
```

or, equivalently, pack it directly with `packmol(system; ionic_concentration=0.15u"mol/L", margin=20.0, cubic=true, output="system.pdb")`.

Enough cations or anions (of a single sign) are added on top of the bulk
salt to exactly neutralize `system.solute_charge`; a `SolutionBoxUWI` with
default (monovalent) ions can neutralize any integer solute charge. An
`ArgumentError` is raised if the requested charge cannot be exactly
neutralized with the ions' charges (for instance, a solute charge of `-3`
with only doubly-charged cations available). `margin`, `box_sides`, and
`cubic` behave as in `SolutionBoxUS`.
