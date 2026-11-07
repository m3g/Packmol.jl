# packing

The packing: the `packmol` functions of the API, and the packing loop.

| File | Contents |
|---|---|
| `packmol.jl` | `packmol(input_file)`, `packmol()` (file dialog) and `packmol(system)`. Their keyword arguments (the global options of the input files and the options of the run) are listed once, in `_packmol_global_options` and `_packmol_run_options`, and the methods are generated from these lists, so that the options are explicit in the signatures (completed in the REPL); the `packmol` methods of the recipes forward the same options. |
| `packing_loop.jl` | `_packmol`, the packing loop: the initial placement (`../initial_approximation/`), then loops of optimization with SPGBox until the convergence criterion is satisfied (no pair of atoms closer than the tolerance by more than `tolerance_precision`, constraints satisfied to `constraint_precision`). The scale of the radii decays from `radscale` to 1 in `radscale_decay_loops` loops; after the first loop, the bad molecules are moved (`movebad!`) when a loop improves the function little. Also the helpers of the loop: the schedule of the radii, the cutoff of the cell lists, the cell of non-periodic systems, and the objective function of the scaled rotations. |
| `active_set.jl` | The active subsystem: each loop optimizes the molecules that fail the convergence criterion and those near them (within `active_radius`), with the others frozen (`_active_subsystem`, `_optimize_active_subsystem!`). |
| `report.jl` | The report of the loop in the terminal (progress bar, overlaps by pairs of structure types), and the `<output>_FORCED` file written when the packing fails. |
