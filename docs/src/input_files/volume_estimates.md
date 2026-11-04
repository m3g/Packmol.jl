# Volume estimates

Two estimates help to choose the numbers of molecules, or the size of the box, of an
input: the volume that the molecules fill at a given density, and the volume of the region
available to the molecules of each structure type, given its constraints and the fixed
molecules. Neither packs the system. From the command line, with the input file:

```
packmol --guess-volume input.inp               # at 1 g/mL
packmol --guess-volume --density 0.8 input.inp
packmol --free-volume input.inp
```

and in Julia, with [`guess_volume`](@ref) and [`free_volume`](@ref), which take an input
file or a [`PackmolSystem`](@ref). Both are also online calculators, in the browser: the
[volume guesser](@ref "Volume guesser") and the [free volume](@ref "Free volume")
calculator, at the top of the Recipes tab.

## The volume of the molecules at a density

`guess_volume` sums the masses of the molecules (the fixed ones included) and divides it
by the density of the system, as the volume guesser of the
[Packmol web page](https://m3g.github.io/packmol): the result is the volume of the box,
and the side of a cube of that volume. For the [water and urea](@ref "Water and urea")
tutorial (1000 water and 400 urea molecules), at 1.09 g/mL:

```@example volumes
using Packmol
guess_volume("../tutorials/files/mixture/mixture.inp"; density = 1.09)
nothing # hide
```

which is the cube of 40 Å of the input. The molecules can also be given as structure files
or molar masses (g/mol), with their numbers:

```@example volumes
guess_volume(["../tutorials/files/mixture/water.pdb" => 1000, 60.06 => 400]; density = 1.09)
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
free_volume("../tutorials/files/solvprotein/solvprotein.inp"; npoints = 500_000)
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
free_volume("../tutorials/files/bilayer/bilayer.inp"; npoints = 500_000)
nothing # hide
```

A structure type whose region is not bounded (without periodic boundary conditions, by an
`inside` box, cube, sphere, cylinder or ellipsoid) has no estimate either. The precision of
the estimate is given by its standard error, which decreases with the square root of the
number of points (`npoints`, 1 million by default).

```@docs
guess_volume
free_volume
```
