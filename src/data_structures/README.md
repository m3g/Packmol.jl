# data_structures

The types that describe a system to be packed.

| File | Contents |
|---|---|
| `atoms_and_molecules.jl` | `AtomData` (the molecule, structure type and radius of each atom) and `MoleculePosition` (the variables of a molecule: center of mass and Euler angles). The vector of `MoleculePosition` is reinterpreted as the flat vector of variables of the solver. |
| `StructureType.jl` | `StructureType`, a type of molecule (structure, number of molecules, constraints, options), `structure_type` (its constructor in the Julia API), the reading of structure files, and the expansion of the conformations of a molecule into consecutive structure types. |
| `PackmolSystem.jl` | `PackmolSystem`, the system (structure types, atoms, variables and options), its constructors and `show`, and the scale of the radii of each structure type at the start of the packing (`_radscale`). |
| `system_options.jl` | The options of a system given as keyword arguments of `PackmolSystem` (and set by the global keywords of an input file): the periodic cell, the tolerance, and the others. |
| `FixedParticleSystem.jl` | Cell lists of the fixed atoms, built once, for the distances between free and fixed atoms in the initial placement and in movebad. |

The input file parser, which builds these types from an input file, is in `../io/`.
