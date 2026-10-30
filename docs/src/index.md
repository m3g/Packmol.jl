# Packmol.jl

Packmol creates initial configurations for molecular dynamics simulations, by packing
molecules in regions of space defined by geometric constraints, such that the atoms of
different molecules are kept at a safe distance from each other.

Packmol.jl is the Julia implementation of [Packmol](http://github.com/m3g/packmol). It
reads the input files of Fortran Packmol and runs them with the same command line, and it
can also be used from Julia code, where systems are defined without input files and can
be built from target densities and concentrations (see [Recipes](recipes/index.md)). It uses
[CellListMap.jl](https://github.com/m3g/CellListMap.jl) for the computation of distances
and [SPGBox.jl](https://github.com/m3g/SPGBox.jl) as the optimizer.

## Documentation

- [Tutorials](tutorials/index.md): complete examples, each with an input file for the
  command line and the same system defined in Julia.
- [Input files](input_files/index.md): the global keywords and the `structure ... end structure`
  blocks.
- [Constraints](constraints/index.md): the geometric regions (`box`, `sphere`, `plane`,
  `cylinder`, `ellipsoid`, ...) in which molecules, or some of their atoms, are packed.
- [Periodic boundary conditions](periodic_boundary_conditions/index.md): orthorhombic,
  triclinic, rhombic dodecahedral and truncated octahedral cells.
- [Defining systems in Julia](julia_api/index.md): building and packing a system with
  `structure_type`, `PackmolSystem` and `packmol`, without an input file.
- [Recipes](recipes/index.md): solvated systems and membranes built from target densities and
  concentrations, instead of numbers of molecules.
- [Concentration units](recipes/concentration_units.md): the conversions between concentration
  units (`cconvert`) used by the recipes.
- [Reference](julia_api/reference.md): the docstrings of the public functions.
