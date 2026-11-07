# function_and_gradient

The objective function of the packing and its gradient: the sum of the penalties of the
pairs of atoms of different molecules closer than the sum of their radii, and of the
violations of the constraints. The pairs are found with cell lists (CellListMap.jl), in
parallel.

| File | Contents |
|---|---|
| `fg_output.jl` | `InteratomicDistanceFG`, the output of the evaluation (function, gradient, minimum distance, function of each molecule, overlaps), and its interface with CellListMap (copy, reset and reduction of the threaded copies). |
| `pair_fg.jl` | The function and gradient of one pair of atoms (`cartesian_fg!`), accumulated per molecule and per atom. |
| `atom_positions.jl` | The positions of the atoms from the variables of the molecules (with the constraints evaluated in the same pass), and their copy into the cell lists, wrapped into the periodic cell. |
| `constraint_fg.jl` | The function and gradient of the constraints of the atoms and of the rotations. |
| `fg.jl` | `fg!`, the objective function given to the solver (with the gradient relative to the variables of the molecules), and `_mark_unconverged!` (the molecules that fail the convergence criterion). |
