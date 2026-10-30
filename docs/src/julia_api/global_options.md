# Global options

Any [`PackmolSystem`](@ref) field can be passed as a keyword to the
`PackmolSystem` constructor — `tolerance`, `seed`, `radscale`,
`tolerance_precision`, `unitcell`, and so on, matching the input file's
[Global keywords](@ref) of the same name. [`packmol`](@ref) itself also takes
a few keywords directly (`nloop`, `maxit`, `iprint`, `parallel`,
...) that don't have an input-file equivalent, since they control this
particular run rather than being system properties.
