# Lipid bilayer

A lipid bilayer with water above and below it. The lipids of each layer are packed in a
slab, and constraints on some of their atoms orient them: the polar heads point to the
water and the hydrophobic tails to the center of the bilayer.

![Lipid bilayer](https://m3g.github.io/packmol/examples/bilayer.jpg)

Files: [`palmitoil.pdb`](https://m3g.github.io/packmol/examples/palmitoil.pdb),
[`water.pdb`](https://m3g.github.io/packmol/examples/water.pdb) and the input file
[`bilayer.inp`](https://m3g.github.io/packmol/examples/bilayer-comment.inp).

In `palmitoil.pdb`, atoms 31 and 32 belong to the polar head of the lipid, and atoms 1
and 2 are the end of its hydrophobic tail.

## Command line

The input file:

```
#
# Lipid double layer with water above and below
#

# Atoms from different molecules will be at least 2.0 Angstroms apart
# in the solution.

tolerance 2.0

# Coordinate files are in the pdb format (the default)

filetype pdb

# The output pdb file

output bilayer.pdb

# Water molecules below the lipids:
# The first three numbers are the minimum x, y, z coordinates of these
# molecules, and the last three are the maximum coordinates. The box
# has sides of 40 Angstroms in the x and y directions, and of
# 10 Angstroms in the z direction.

structure water.pdb
  number 500
  inside box 0. 0. -10. 40. 40. 0.
end structure

# Water molecules above the lipids:
# The same as above, but with the box in a different region of space.

structure water.pdb
  number 500
  inside box 0. 0. 28. 40. 40. 38.
end structure

# First lipid layer: the polar heads point down, to the water below.
# 50 lipids will be put inside a box of side 40 in the x and y
# directions and 14 in the z direction (a little more than the length
# of the lipid). Atoms 31 and 32, of the polar head, are constrained to
# be below the plane z = 2, and atoms 1 and 2, the hydrophobic end of
# the lipid, above the plane z = 12. Therefore, all the lipids will be
# oriented, with their polar heads pointing to the water box below.

structure palmitoil.pdb
  number 50
  inside box 0. 0. 0. 40. 40. 14.
  atoms 31 32
    below plane 0. 0. 1. 2.
  end atoms
  atoms 1 2
    above plane 0. 0. 1. 12.
  end atoms
end structure

# Second lipid layer: the polar heads point up, to the water above,
# with the opposite orientation of the first layer.

structure palmitoil.pdb
  number 50
  inside box 0. 0. 14. 40. 40. 28.
  atoms 1 2
    below plane 0. 0. 1. 16.
  end atoms
  atoms 31 32
    above plane 0. 0. 1. 26
  end atoms
end structure
```

The constraints inside an `atoms ... end atoms` block apply only to the atoms listed,
while those outside it (here, `inside box`) apply to all the atoms of the molecule.
`below plane 0. 0. 1. 2.` is the half-space `0x + 0y + 1z ≤ 2`, that is, `z ≤ 2`.

Run it with:

```
packmol bilayer.inp
```

The packed system is written to `bilayer.pdb`.

## Julia

The `atom_constraints` keyword of [`structure_type`](@ref) is the equivalent of the
`atoms ... end atoms` blocks: a list of `atom indices => constraints` pairs.

```julia
using Packmol

water_below = structure_type("water.pdb";
    number = 500,
    constraints = [InsideBox([0., 0., -10.], [40., 40., 0.])],
)
water_above = structure_type("water.pdb";
    number = 500,
    constraints = [InsideBox([0., 0., 28.], [40., 40., 38.])],
)
# Polar heads (atoms 31 and 32) down, tails (atoms 1 and 2) up
lower_layer = structure_type("palmitoil.pdb";
    number = 50,
    constraints = [InsideBox([0., 0., 0.], [40., 40., 14.])],
    atom_constraints = [
        [31, 32] => BelowPlane([0., 0., 1.], 2.0),
        [1, 2] => AbovePlane([0., 0., 1.], 12.0),
    ],
)
# Tails down, polar heads up
upper_layer = structure_type("palmitoil.pdb";
    number = 50,
    constraints = [InsideBox([0., 0., 14.], [40., 40., 28.])],
    atom_constraints = [
        [1, 2] => BelowPlane([0., 0., 1.], 16.0),
        [31, 32] => AbovePlane([0., 0., 1.], 26.0),
    ],
)

system = PackmolSystem([water_below, water_above, lower_layer, upper_layer];
    output = "bilayer.pdb",
    tolerance = 2.0,
)
packmol(system)
```

`BelowPlane(normal, d)` and `AbovePlane(normal, d)` take the parameters of
`below plane` and `above plane`.

!!! tip
    The `Membrane` recipe builds bilayers like this one from a few parameters
    (lipid composition, area per lipid, solvent thickness), with periodic boundary
    conditions. See [Recipes](../experimental/recipes.md).
