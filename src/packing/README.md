# packing

The packing: the `packmol` functions of the API, and the packing loop.

| File | Contents |
|---|---|
| `packmol.jl` | `packmol(input_file)`, `packmol()` (file dialog) and `packmol(system)`, with their keyword arguments (`nloop`, `maxit`, `nthreads`, ...). |
| `packing_loop.jl` | `_packmol`, the packing loop: the initial placement (`../initial_approximation/`), then loops of optimization with SPGBox until the convergence criterion is satisfied (no pair of atoms closer than the tolerance by more than `tolerance_precision`, constraints satisfied to `constraint_precision`). The scale of the radii decays from `radscale` to 1 in `radscale_decay_loops` loops; after the first loop, the bad molecules are moved (`movebad!`) when a loop improves the function little. Also the helpers of the loop: the schedule of the radii, the cutoff of the cell lists, the cell of non-periodic systems, and the objective function of the scaled rotations. |
| `active_set.jl` | The active subsystem: each loop optimizes the molecules that fail the convergence criterion and those near them (within `active_radius`), with the others frozen (`_active_subsystem`, `_optimize_active_subsystem!`). |
| `report.jl` | The report of the loop in the terminal (progress bar, overlaps by pairs of structure types), and the `<output>_FORCED` file written when the packing fails. |
