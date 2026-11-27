```@meta
CollapsedDocStrings = true
```
# Volume estimates

Three estimates help to choose the numbers of molecules, or the size of the box, of an
input: the volume that the molecules fill at a given density, the volume of the region
available to the molecules of each structure type, given its constraints and the fixed
molecules, and the volume of the molecules themselves. None packs the system. From the
command line, with the input file:

```
packmol --guess-volume input.inp               # at 1 g/mL
packmol --guess-volume --density 0.8 input.inp
packmol --free-volume input.inp
packmol --molecular-volume input.inp
packmol --molecular-volume --resolution 0.2 input.inp
```

and in Julia, with [`guess_volume`](@ref), [`free_volume`](@ref) and
[`molecular_volume`](@ref), which take an input file or a [`PackmolSystem`](@ref). Both are also online calculators, in the browser: the
[volume guesser](@ref "Volume guesser") and the [free volume](@ref "Free volume")
calculator, at the top of the Recipes tab.

## The volume of the molecules at a density

`guess_volume` sums the masses of the molecules (the fixed ones included) and divides it
by the density of the system, as the volume guesser of the
[Packmol web page](https://m3g.github.io/packmol): the result is the volume of the box,
and the side of a cube of that volume. For the [water and urea](@ref "Water and urea")
tutorial (1000 water and 400 urea molecules), at 1.09 g/mL:

```@setup volumes
# The input files of the tutorials used below, each in a directory with the structure
# files of the tutorials (the packmol_tutorials set of MolSimTestFiles)
using Packmol
for t in ("mixture", "solvprotein", "bilayer")
    Main.copy_tutorial(t, t)
end
```

```@example volumes
using Packmol
guess_volume("mixture/mixture.inp"; density = 1.09)
nothing # hide
```

which is the cube of 40 Å of the input. The molecules can also be given as structure files
or molar masses (g/mol), with their numbers:

```@example volumes
water = "mixture/water.pdb"
guess_volume([water => 1000, 60.06 => 400]; density = 1.09)
nothing # hide
```

## The volume available to each structure type

`free_volume` estimates, by Monte Carlo, the volume of the region of each structure type:
the fraction of random points of a box that contains the region that satisfy its
constraints and are farther than the tolerance from the atoms of the fixed molecules. With
it comes the density at which the molecules of the structure type fill their region: much
above the density of the liquid, there are too many molecules for the region (and the
packing will be slow, or fail); much below, too few.

For the [solvated protein](@ref "Solvated protein"), the water, sodium and chloride share
the rhombic dodecahedron, without the volume of the protein:

```@example volumes
input = "solvprotein/solvprotein.inp"
free_volume(input; npoints = 500_000)
nothing # hide
```

The water fills its region at 1 g/mL. The densities of structure types that share a region
add up: the water and the urea of the mixture above, in the same box, are at 0.47 and 0.62
g/mL, 1.09 g/mL together.

The region is that of the constraints that apply to the whole molecule. The constraints of
some atoms only (in `atoms` blocks) are ignored, so that the volume is an upper bound, and a
structure type with constraints of some atoms only has no estimate: in the [lipid
bilayer](@ref "Lipid bilayer"), the heads and the ends of the chains of the lipids are above
and below planes, which do not bound a region of the molecules. The water around them has
its own region:

```@example volumes
free_volume("bilayer/bilayer.inp"; npoints = 500_000)
nothing # hide
```

A structure type whose region is not bounded (without periodic boundary conditions, by an
`inside` box, cube, sphere, cylinder or ellipsoid) has no estimate either. The precision of
the estimate is given by its standard error, which decreases with the square root of the
number of points (`npoints`, 1 million by default).

## The volume of the molecules

`molecular_volume` computes the volume of each molecule from its atoms, as
[MoloVol](https://molovol.com) does (Maglic & Lavendomme, J. Appl. Cryst. 55, 1033
(2022)): the volume from which a spherical probe (of 1.4 Å, the radius of a water
molecule) is excluded by the van der Waals spheres of the atoms, on a grid. It is the van
der Waals volume of the atoms and the gaps between them where the probe does not fit;
cavities where the probe fits are not part of it.

It is the volume that a large solute takes from the solvent: for a protein, close to the
volume given by its partial specific volume in water (about 0.73 mL/g), while its mass at
the density of water gives about 1.4 times that. The solvation [recipes](@ref "Recipes")
subtract it from the volume of the box. For the protein of the [solvated
protein](@ref "Solvated protein") tutorial:

```@example volumes
molecular_volume(["solvprotein/protein.pdb" => 1])
nothing # hide
```

The grid has a spacing of 0.4 Å by default (`resolution`), at which the volume is about
2% larger than the converged one, in a fraction of a second for thousands of atoms. A
smaller spacing gives a more precise volume, at a cost that grows as its inverse cube.
The radii are the van der Waals radii of the elements, so that the volume is that of
all-atom structures. With `probe = 0`, it is the van der Waals volume.

For small molecules, the molecular volume is smaller than their volume in the liquid (by
20 to 40%: it does not include the empty space between the molecules of a liquid), which
the density of the liquid gives. For the solvent, use [`guess_volume`](@ref).

```@docs
guess_volume
free_volume
molecular_volume
```
