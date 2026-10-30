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

## Example

830 water molecules (about 1 g/mL) in a triclinic cell, with sides of 30 Å and angles of 70°, 80° and 90°:

```
unitcell 30. 30. 30. 70. 80. 90.

structure water.pdb
    number 830
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl):

```@setup triclinic
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example triclinic
using Packmol, PDBTools
water = structure_type("water.pdb"; number = 830)
system = PackmolSystem([water]; tolerance = 2.0, output = "triclinic.pdb",
    unitcell = UnitCell(30., 30., 30., 70., 80., 90.))
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system); unitcell = UnitCell(system.unitcell), unitcell_origin = -sum(eachcol(system.unitcell)) / 2)
```
