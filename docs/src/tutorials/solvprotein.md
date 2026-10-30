# Solvated protein

A protein of 4246 atoms in a sphere of water molecules and sodium and chloride ions. The
protein is fixed at the center of the sphere, and the solvent fills the space around it.
This system was used in [Martínez L., et al. J. Med. Chem. 49, 23–26, 2006](https://doi.org/10.1021/jm050805n).

![Solvated protein](https://m3g.github.io/packmol/examples/solvprotein.jpg)

Files: [`protein.pdb`](https://m3g.github.io/packmol/examples/protein.pdb),
[`water.pdb`](https://m3g.github.io/packmol/examples/water.pdb),
[`CLA.pdb`](https://m3g.github.io/packmol/examples/CLA.pdb),
[`SOD.pdb`](https://m3g.github.io/packmol/examples/SOD.pdb) and the input file
[`solvprotein.inp`](https://m3g.github.io/packmol/examples/solvprotein-comment.inp).

## Command line

The input file:

```
#
# Protein solvated by water and ions
#

# All atoms from different molecules will be at least 2.0 Angstroms apart
# in the solution

tolerance 2.0

# The type of the files will be pdb

filetype pdb

# The name of the output file

output solvprotein.pdb

# The protein will be fixed with its center of mass at the center of the
# sphere, and no rotation (the first three zeros are the position of the
# center of mass and the last three are the rotation angles, in radians,
# relative to the orientation in the input file).

structure protein.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  centerofmass
end structure

# Water molecules will be put inside a sphere that contains the protein:
# 16500 water molecules in a sphere centered at the origin, with a
# radius of 50 Angstroms.

structure water.pdb
  number 16500
  inside sphere 0. 0. 0. 50.
end structure

# The ions will be put in the same region as the water molecules, so
# they will also be solvating the protein.

structure CLA.pdb
  number 20
  inside sphere 0. 0. 0. 50.
end structure

structure SOD.pdb
  number 30
  inside sphere 0. 0. 0. 50.
end structure
```

No constraint keeps the solvent away from the protein: the protein is fixed, and the
packing itself keeps the other molecules at least `tolerance` away from its atoms.

Run it with:

```
packmol solvprotein.inp
```

The packed system is written to `solvprotein.pdb`.

## Julia

```julia
using Packmol

protein = structure_type("protein.pdb";
    number = 1,
    fixed = ([0., 0., 0.], [0., 0., 0.]),
    center = :mass,
)
sphere = InsideSphere([0., 0., 0.], 50.0)
water = structure_type("water.pdb"; number = 16500, constraints = [sphere])
chloride = structure_type("CLA.pdb"; number = 20, constraints = [sphere])
sodium = structure_type("SOD.pdb"; number = 30, constraints = [sphere])

system = PackmolSystem([protein, water, chloride, sodium];
    output = "solvprotein.pdb",
    tolerance = 2.0,
)
packmol(system)
```

A constraint is a regular Julia value, so the same `sphere` can be given to several
structure types.

!!! tip
    The `SolutionBoxUWI` recipe solvates a protein with water and ions from a
    target ionic concentration (neutralizing its charge), in a periodic box. See
    [Recipes](../recipes/index.md).
