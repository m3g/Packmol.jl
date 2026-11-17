```@meta
CollapsedDocStrings = true
```
# Recipes

A recipe is a higher-level, parameter-driven system setup: instead of
specifying molecule counts directly (as a `structure ... end structure`
block does), a recipe takes target densities or concentrations and works out
the molecule counts and box size itself. Two online calculators, the
[volume guesser](@ref "Volume guesser") and the [free volume](@ref "Free volume")
calculator, help to choose them for any input. Recipes currently cover common
solvation setups:

- [Solute-Solvent system](@ref): a solute solvated by a single solvent.
- [Solute-Solvent-Cossolvent system](@ref): a solute solvated by a mixture of
  two solvents (e.g. water/ethanol), at a target cossolvent concentration.
- [Solute-Water-Ions system](@ref): a solute solvated by water and a
  background salt (e.g. NaCl) at a target ionic concentration, with
  automatic charge neutralization of the solute.
- [Membrane system](@ref): a lipid bilayer or monolayer of one or more
  lipids, solvated above (and, for a bilayer, also below) the membrane. Seven
  common phospholipids are built in: only the composition and the size are needed.
- [Membrane protein system](@ref): a membrane protein (oriented as in the OPM
  database) in a lipid bilayer of a given composition, solvated by water and
  ions.

Vesicles and nanotubes are planned as future recipes.

Concentrations and densities can be given as plain numbers (assumed to be in
the units noted in each function's docstring) or as [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/)
quantities (`55.5u"mol/L"`, `1.0u"g/mL"`, ...). See
[Concentration Unit Conversion](@ref) for the underlying `cconvert` machinery
used to convert between them.

## [Packing a recipe](@id recipes_packing)

A recipe is a data structure describing the components of the system (structure files,
molar masses, densities), and `packmol` builds and packs it directly, without an input
file:

```julia
using Packmol
system = SolutionBoxUS(...)
packmol(system; margin = 20.0, output = "system.pdb")
```

This builds a `PackmolSystem` in memory (the structure of the [Julia
API](../julia_api/index.md)) and packs it with `packmol(::PackmolSystem)`: `output` is the
only file written, and the packed system is returned. The keywords of each recipe
(`margin`/`box_sides`/`pbc`, and `concentration`/`ionic_concentration` where relevant) are
documented in its page; any other keyword (`nloop`, `maxit`, `seed`, ...) is passed to the
packing engine (see `packmol(::PackmolSystem)`).

### The periodic cell

In the solvation recipes, `pbc` selects the shape of the periodic cell, sized by `margin`
or `box_sides`: `:cubic` (the default), `:orthorhombic`, `:dodecahedral` (a rhombic
dodecahedron) or `:octahedral` (a truncated octahedron). See [Cell shapes](@ref
cell_shapes) for how each is sized, and the pages of the shapes. The numbers of molecules
follow from the volume of the cell, whatever its shape. `Membrane` and `MembraneProtein`
are always orthorhombic.

## [Writing a Packmol input file](@id recipes_input)

`write_packmol_input` writes the input file of Packmol for a recipe, to be inspected,
edited or run later, with the `packmol` executable or in Julia:

```julia
using Packmol
write_packmol_input(system;
    margin = 20.0,
    input = "box.inp",
    output = "system.pdb",
)
packmol("box.inp")
```

The input file can also be run with the Fortran Packmol, wrapped by
[`run_packmol_legacy`](@ref):

```julia
using Packmol, Packmol_jll
run_packmol_legacy("box.inp")
```
