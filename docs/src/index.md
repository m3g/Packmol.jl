# Packmol.jl

Packmol creates initial configurations for molecular dynamics simulations, by packing
molecules in regions of space defined by geometric constraints, such that the atoms of
different molecules are kept at a safe distance from each other.

Packmol.jl is the new implementation of [Packmol](https://doi.org/10.1002/jcc.21224). It
reads the input files of the legacy (Fortran) Packmol and runs them with the same command
line, and it can also be used from Julia code, where systems are defined without input
files and can be built from target densities and concentrations (see
[Recipes](recipes/index.md)). It uses
[CellListMap.jl](https://github.com/m3g/CellListMap.jl) for the computation of distances
and [SPGBox.jl](https://github.com/m3g/SPGBox.jl) as the optimizer.

## New in Packmol.jl

Relative to the legacy (Fortran) Packmol:

- **Much faster, which makes new systems possible.** Packmol.jl is 2 to 16 times faster in
  a single thread, and it is multithreaded. Systems that were out of reach become routine:
  nanoliposomes of 350 and 600 thousand atoms are packed in 16 and 31 seconds, a whole
  [Zika virion](tutorials/zika.md) (2.2 million particles) in under a minute, an
  [exosome](tutorials/exosome.md) of 13.6 million atoms in about 26 minutes, and the
  [minimal cell](tutorials/cell.md) of Bentopy, with about 40 million beads, in about 11
  minutes (14 threads).
- **Periodic boundary conditions** in orthorhombic, triclinic, rhombic dodecahedral and
  truncated octahedral cells.
- **New constraints**: curved surfaces (Gaussian, waves, exponentials, and their radial
  versions), surfaces given by points (`within surface`, `inside surface`), and
  `constrain_rotation` for the angle between an axis and a vector of the molecule (to
  build oriented phases, such as liquid crystals).
- **Crowded systems**: [`voxel_place`](@ref), a direct port to Julia of the placement
  algorithm of [Bentopy](https://github.com/marrink-lab/bentopy), for the initial
  placement of very crowded systems, such as the cytosol of a cell.
- **mmCIF** files, for structures, outputs and restarts, and restarts from the output of a
  previous run.
- **Several conformations** of a molecule, from a structure file with several models.
- **Volume estimates**, to compute the numbers of molecules (`--guess-volume`,
  `--free-volume`), and `check_initial_placement`, to check the constraints before packing.
- **A Julia API**: systems defined without input files, **recipes** for solvated systems
  and membranes built from target densities and concentrations, user-defined constraints
  ([`CustomConstraint`](@ref)), and the packed atoms in memory, for analysis with
  [PDBTools.jl](https://github.com/m3g/PDBTools.jl). It is also the best way to build
  systems with an [AI agent](tutorials/ai_agent.md).
- **A standalone executable** for Linux, macOS and Windows, installed with `pip`/`uvx` or
  downloaded (see [Installation](installation.md)).

The input files of the legacy Packmol run unchanged, with few exceptions: see
[Migrating from legacy Packmol](migrating.md).

## Documentation

- [Migrating from legacy Packmol](migrating.md): what changed from the Fortran Packmol, what
  is no longer supported, and what is new.
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

## References

If Packmol contributed to your work, please cite the following references:

- L. Martínez, R. Andrade, E. G. Birgin, J. M. Martínez. Packmol: A package for building
  initial configurations for molecular dynamics simulations. *Journal of Computational
  Chemistry*, 30(13):2157-2164, 2009. [doi:10.1002/jcc.21224](https://doi.org/10.1002/jcc.21224)

- J. M. Martínez and L. Martínez. Packing optimization for automated generation of complex
  system's initial configurations for molecular dynamics and docking. *Journal of
  Computational Chemistry*, 24(7):819-825, 2003. [doi:10.1002/jcc.10216](https://doi.org/10.1002/jcc.10216)

## See also

- [ComplexMixtures.jl](https://m3g.github.io/ComplexMixtures.jl/stable/): solvation and
  interactions in complex solutions (minimum-distance distribution functions, Kirkwood-Buff
  integrals), from simulations.
- [CellListMap.jl](https://m3g.github.io/CellListMap.jl/stable/): fast computation of
  pairwise interactions and distances with cell lists, on which Packmol.jl is built.
- [PDBTools.jl](https://m3g.github.io/PDBTools.jl/stable/): reading, writing, selecting and
  visualizing molecular structures (PDB and mmCIF).
- [m3g.github.io](https://m3g.github.io/): the other packages and the research of our group.

Contact: [lmartine@unicamp.br](mailto:lmartine@unicamp.br)
