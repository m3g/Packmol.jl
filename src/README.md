# Source code of Packmol.jl

`Packmol.jl` defines the module: its dependencies, the exported names, and the list of
included files, in an order in which every type is defined before it is used in a method
signature. Each directory has a README of its own.

| Directory | Contents |
|---|---|
| [`utils/`](utils/README.md) | Terminal output, progress bar, task synchronization (all compatible with `juliac --trim`), and the file dialog |
| [`constraints/`](constraints/README.md) | The geometric constraints (boxes, spheres, planes, curved surfaces, custom functions, rotations): penalties, gradients and parsing |
| [`data_structures/`](data_structures/README.md) | `StructureType`, `PackmolSystem` and the data of atoms and molecules |
| [`periodic_cells/`](periodic_cells/README.md) | Periodic cells (orthorhombic, triclinic, rhombic dodecahedron, truncated octahedron), and wrapping into them |
| [`io/`](io/README.md) | The input file parser, the output files, and restarts |
| [`rigid_body/`](rigid_body/README.md) | Rotations and translations of rigid molecules, the chain rule of the gradient, and the scaling of the rotations |
| [`function_and_gradient/`](function_and_gradient/README.md) | The objective function of the packing and its gradient |
| [`initial_approximation/`](initial_approximation/README.md) | The initial placement of the molecules, and the relocation of bad molecules (movebad) |
| [`packing/`](packing/README.md) | The `packmol` functions and the packing loop, with the active subsystem |
| [`tools/`](tools/README.md) | Tools for setting up systems: volumes, shapes of molecules, lipids and cross-sections |
| [`Recipes/`](Recipes/README.md) | Systems set up from densities and concentrations (solvation boxes, membranes) |

Files at this level:

- `Packmol.jl`: the module.
- `CLI.jl`: the command line interface (the `packmol` app and the executable built with
  `juliac`, see `juliac/`).
- `packmol_runner.jl`: `run_packmol_legacy`, the runner of the Fortran Packmol (its methods
  are in the `PackmolJLLExt` extension, `ext/`).
- `precompile.jl`: the precompilation workload (skipped when the executable is built).

The tests (`test/`) mirror this layout: the tests of `src/<dir>/<file>.jl` are in
`test/<dir>/<file>.jl`.

## Conventions

- Each file starts with a comment saying what it contains.
- Dimensions and precision are type parameters (`D`, `T`); 3D is the main target.
- Coordinates are in Å.
- Names starting with `_` are internal. Keyword arguments passed through are `kwargs...`.
- Errors are `ArgumentError`s with messages padded as `"""\n ... \n\n"""`.
- The code reachable from the command line interface must compile with
  `juliac --trim=safe`: see the notes in `CLAUDE.md` (no printing to `stdout` with
  `println`, no `@sync`, concrete types in closures, ...).
