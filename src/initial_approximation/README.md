# initial_approximation

The initial placement of the molecules (`set_initial_approximation!`, following the
steps of Fortran Packmol), and the relocation of bad molecules during the packing
(`movebad!`):

1. random positions, in the region of the constraints or in the periodic cell
   (`initialize_molecules.jl`);
2. a fit of the positions to the constraints alone (`adjust_constraints.jl`), from which
   the bounds of the centers of mass of each structure type are computed;
3. a new placement, the best of many random trials for each molecule, within those bounds
   (`reinitialize_with_bounds.jl`);
4. a fit to the constraints and to the distances to the fixed molecules, with the
   molecules that do not satisfy them moved;
5. the molecules still far outside the bounds are clamped to them (`bounds.jl`).

With `check_initial_placement`, the result of these steps is written, and the packing
stops.

| File | Contents |
|---|---|
| `initial_approximation.jl` | `set_initial_approximation!`, the sequence of steps. |
| `initialize_molecules.jl` | Step 1. |
| `adjust_constraints.jl` | The fits to the constraints (and to the distances to the fixed atoms). |
| `reinitialize_with_bounds.jl` | Step 3. |
| `randomize_molecule.jl` | A random position and orientation of a molecule, within bounds or in the periodic cell. |
| `molecule_badness.jl` | The badness of each molecule (constraint penalty and overlap with the fixed molecules), to select the best trials. |
| `movebad.jl` | `movebad!`: a fraction of the molecules that fail the convergence criterion moved to random positions, in the packing loops (after the first) that improve the function little. |
| `bounds.jl` | Bounds of the atom positions and of the centers of mass of each structure type, and the clamping of molecules to them. |
| `constraint_checks.jl` | The constraints of single molecules: whether they are satisfied, their penalty, and the function and gradient of the constraints alone. |
| `fixed_molecules.jl` | The atoms of the fixed molecules, and the overlaps of the free molecules with them. |
| `buffers.jl` | Buffers and index arrays shared by these functions. |
