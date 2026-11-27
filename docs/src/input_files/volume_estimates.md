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
packmol --guess-volume input.inp
packmol --free-volume input.inp
packmol --molecular-volume input.inp
```

and in Julia, with [`guess_volume`](@ref), [`free_volume`](@ref) and
[`molecular_volume`](@ref), which print the same summaries. The first two are also
online calculators, in the browser: the [volume guesser](@ref "Volume guesser") and the
[free volume](@ref "Free volume") calculator, at the top of the Recipes tab.

```@setup volumes
# The input files of the tutorials used below, each in a directory with the structure
# files of the tutorials (the packmol_tutorials set of MolSimTestFiles)
using Packmol
for t in ("mixture", "solvprotein", "bilayer")
    Main.copy_tutorial(t, t)
end
```

## The volume of the molecules at a density

The volume of the box, and the side of a cube of that volume: the mass of the free
molecules divided by the density, as the volume guesser of the
[Packmol web page](https://m3g.github.io/packmol), plus the volume of the fixed molecules,
their [molecular volume](@ref "The volume of the molecules") (as the solvation recipes
count the volume of the solute).

From the command line, `--density` sets the density (g/mL, 1.0 by default):

```
packmol --guess-volume --density 1.09 input.inp
```

In Julia, `guess_volume` takes an input file, a [`PackmolSystem`](@ref), or structure
files or molar masses (g/mol) with their numbers of molecules, and the `density` keyword.
For the [water and urea](@ref "Water and urea") tutorial (1000 water and 400 urea
molecules), at 1.09 g/mL:

```@example volumes
using Packmol
guess_volume("mixture/mixture.inp"; density = 1.09)
nothing # hide
```

which is the cube of 40 Å of the input. The same, with the molecules given directly (all
free):

```@example volumes
water = "mixture/water.pdb"
guess_volume([water => 1000, 60.06 => 400]; density = 1.09)
nothing # hide
```

For the [solvated protein](@ref "Solvated protein"), the protein is fixed, and counts
with its molecular volume:

```@example volumes
guess_volume("solvprotein/solvprotein.inp")
nothing # hide
```

## The volume available to each structure type

The volume of the region of each structure type, by Monte Carlo: the fraction of random
points of a box that contains the region that satisfy its constraints and are farther
than the tolerance from the atoms of the fixed molecules. With it comes the density at
which the molecules of the structure type fill their region: much above the density of
the liquid, there are too many molecules for the region (and the packing will be slow, or
fail); much below, too few.

The region is that of the constraints that apply to the whole molecule. The constraints of
some atoms only (in `atoms` blocks) are ignored, so that the volume is an upper bound, and a
structure type with constraints of some atoms only has no estimate. A structure type whose
region is not bounded (without periodic boundary conditions, by an `inside` box, cube,
sphere, cylinder or ellipsoid) has no estimate either. The precision of the estimate is
given by its standard error, which decreases with the square root of the number of points.

From the command line, `--npoints` sets the number of points (1 000 000 by default):

```
packmol --free-volume --npoints 500000 input.inp
```

In Julia, `free_volume` takes an input file or a [`PackmolSystem`](@ref), and the `npoints`
keyword. For the [solvated protein](@ref "Solvated protein"), the water, sodium and
chloride share the rhombic dodecahedron, without the volume of the protein:

```@example volumes
input = "solvprotein/solvprotein.inp"
free_volume(input; npoints = 500_000)
nothing # hide
```

The water fills its region at 1 g/mL. The densities of structure types that share a region
add up: the water and the urea of the mixture above, in the same box, are at 0.47 and 0.62
g/mL, 1.09 g/mL together. In the [lipid bilayer](@ref "Lipid bilayer"), the heads and the
ends of the chains of the lipids are above and below planes (constraints of some atoms),
which do not bound a region of the molecules; the water around them has its own region:

```@example volumes
free_volume("bilayer/bilayer.inp"; npoints = 500_000)
nothing # hide
```

## The volume of the molecules

The volume of each molecule, as the space that a spherical probe cannot enter, computed on
a grid as by [MoloVol](https://molovol.com). With the default probe (1.4 Å, a water
molecule), it is the volume that a solute takes from the solvent: for a protein, close to
its volume in water, while its mass at 1 g/mL gives about 1.4 times that. It is the
volume of the fixed molecules in the [volume at a density](@ref "The volume of the molecules at a density"),
and that of the solute in the solvation [recipes](../recipes/index.md). For small
molecules, as the solvent, the volume in the liquid is larger (it includes the empty space
between the molecules): the volume at a density is the one for them.

The result lists, for each structure of the input file, the volume of one molecule and of
all its molecules, and the total. From the command line:

```
packmol --molecular-volume input.inp
```

`--probe` sets the radius of the probe (Å, 1.4 by default; `0` gives the van der Waals
volume), and `--resolution` the spacing of the grid (Å, 0.4 by default, at which the
volume is within about 1.5% of its converged value).

In Julia, `molecular_volume` takes an input file, a [`PackmolSystem`](@ref), or structure
files with their numbers of molecules, and the `probe` and `resolution` keywords. For the
[solvated protein](@ref "Solvated protein"):

```@example volumes
molecular_volume("solvprotein/solvprotein.inp")
nothing # hide
```

It also gives the volume of a single molecule, from its atoms:

```@example volumes
using PDBTools
molecular_volume(read_pdb("solvprotein/thyroid_receptor.pdb"))
```

```@docs
guess_volume
free_volume
molecular_volume
```
