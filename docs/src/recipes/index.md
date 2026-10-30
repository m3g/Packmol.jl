```@meta
CollapsedDocStrings = true
```
# Recipes

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

## How to use it

```julia-repl
julia> using Packmol
```

`SolutionBoxUS`, `SolutionBoxUSC`, `SolutionBoxUWI`, `Membrane`,
`write_packmol_input`, `convert_concentration`, `convert_density_table!` and
`get_density_table` are exported directly from `Packmol`.
