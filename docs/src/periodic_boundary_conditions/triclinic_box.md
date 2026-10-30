# Triclinic box

```
unitcell a b c alpha beta gamma
```

The fully general periodic cell, CRYST1-style: side lengths `a`, `b`, `c`
and angles `alpha`, `beta`, `gamma` (in degrees, between the corresponding
pairs of sides), built with the PDB convention (`a` along x, `b` in the
xy-plane). Always centered at the origin from the input file; the Julia API
accepts any `unitcell_center` alongside an explicit 3×3 `unitcell` matrix
built however the caller likes (`dodecahedral_unitcell`/`octahedral_unitcell`
below build exactly this kind of matrix from a single size parameter, as a
convenience).

An orthorhombic `pbc` is just the special case `alpha = beta = gamma = 90°`
— everything below (minimum-image distances, implicit confinement, CRYST1
output) is expressed in terms of the general triclinic cell and reduces to
the simpler orthorhombic behavior automatically when the cell happens to be
diagonal.
