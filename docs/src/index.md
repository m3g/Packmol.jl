# Packmol.jl

Packmol creates initial configurations for molecular dynamics simulations, by packing
molecules in regions of space defined by geometric constraints, such that the atoms of
different molecules are kept at a safe distance from each other.

Packmol.jl is the new implementation of [Packmol](http://github.com/m3g/packmol). It
reads the input files of the legacy (Fortran) Packmol and runs them with the same command
line, and it can also be used from Julia code, where systems are defined without input
files and can be built from target densities and concentrations (see
[Recipes](recipes/index.md)). It uses
[CellListMap.jl](https://github.com/m3g/CellListMap.jl) for the computation of distances
and [SPGBox.jl](https://github.com/m3g/SPGBox.jl) as the optimizer.

## New in Packmol.jl

Relative to the legacy (Fortran) Packmol:

- **Much faster**: 3 to 16 times faster in a single thread, and multithreaded, up to about
  40 times faster on a multicore computer (see the table below). Large systems are packed in
  minutes: nanoliposomes of 350 and 600 thousand atoms in 3 and 10 minutes, and a whole
  Zika virion (2.2 million particles) in 7 minutes, with 14 threads.
- **Periodic boundary conditions** in orthorhombic, triclinic, rhombic dodecahedral and
  truncated octahedral cells.
- **New constraints**: curved surfaces (Gaussian, waves, exponentials), and
  `constrain_rotation` for the angle between an axis and a vector of the molecule (to
  build oriented phases, such as liquid crystals).
- **mmCIF** files, for structures, outputs and restarts, and restarts from the output of a
  previous run.
- **Several conformations** of a molecule, from a structure file with several models.
- **Volume estimates**, to compute the numbers of molecules (`--guess-volume`,
  `--free-volume`), and `check_initial_placement`, to check the constraints before packing.
- **A Julia API**: systems defined without input files, **recipes** for solvated systems
  and membranes built from target densities and concentrations, user-defined constraints,
  and the packed atoms in memory, for analysis with
  [PDBTools.jl](https://github.com/m3g/PDBTools.jl).

Packing times of the examples of the legacy Packmol (seconds; the executable of
Packmol.jl, and of Packmol 21.1.4, on an Intel Core i9-14900KF):

| Input | Legacy Packmol | Packmol.jl, 1 thread | Packmol.jl, 14 threads |
|:--|--:|--:|--:|
| water and urea (`mixture`) | 1.4 | 0.14 | 0.16 |
| solvated protein (`solvprotein`) | 13.1 | 1.1 | 0.48 |
| lipid bilayer (`bilayer`) | 26.1 | 1.7 | 0.73 |
| spherical vesicle (`spherical`) | 36.3 | 10.2 | 0.92 |

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

## About the logo

The logo is an image of the [Zika virion](tutorials/zika.md) built with Packmol.jl: a
complete virion of the Zika virus, in a coarse-grained model, with its protein shell, a
membrane of two layers of lipids beneath it, and water and ions inside and around it, in
a rhombic dodecahedron periodic cell. The system has 2.2 million particles, and it is
packed in a few minutes.
