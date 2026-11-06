# Global options

Any [`PackmolSystem`](@ref) field can be passed as a keyword to the
`PackmolSystem` constructor — `tolerance`, `seed`, `radscale`,
`tolerance_precision`, `unitcell`, and so on, matching the input file's
[Global keywords](@ref) of the same name.

Every global keyword of the input files can also be given to [`packmol`](@ref), with the
name of the keyword, and replaces the value of the system (from its input file or its
constructor) before the packing:

```julia
system = read_packmol_input("mixture.inp")
packmol(system; tolerance = 2.5, maxit = 50, seed = 7, output = "mixture_tol2.5.pdb")
```

Options take the values of the input keywords (`periodic_boundary_style = :octahedral`,
`filetype = "mmcif"`), flags `true` or `false` (`add_box_sides = true`), and `pbc` the
values of its input keyword (`pbc = [40, 40, 40]`, `pbc = (:octahedral, 0, 0, 0, 80)`).
[`packmol`](@ref) also takes options of the run that have no input-file equivalent
(`iprint`, `parallel`, `active_set`, ...).
