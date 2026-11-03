# Lipid bilayer

A POPC bilayer in water, with periodic boundary conditions. The lipids of each leaflet
are oriented by constraints on two of their atoms: the polar heads point to the water,
and the ends of the hydrophobic chains to the middle of the bilayer.

![Lipid bilayer](images/bilayer.webp)

*The bilayer seen from the side, with the edges of the periodic cell.*

Files: [`popc.pdb`](files/bilayer/popc.pdb), [`water.pdb`](files/bilayer/water.pdb) and
the input file [`bilayer.inp`](files/bilayer/bilayer.inp).

## The lipid

How well lipids pack depends a lot on their conformation: extended, compact
conformations, with the two chains side by side, pack at the densities of real
membranes, while others do not pack at all. `popc.pdb` is one of the conformations of
the `packmol_lipids` set of
[MolSimTestFiles](https://github.com/m3g/MolSimTestFiles), which has conformations
chosen to pack well for POPC, POPE, POPS, POPG, DPPC, DOPC and DMPC, and describes how
they were obtained (download the whole set:
[packmol_lipids.tar.gz](https://github.com/m3g/MolSimTestFiles/releases/download/packmol_lipids-v1/packmol_lipids.tar.gz)).

In `popc.pdb`, atom 1 is the nitrogen of the choline (the polar head), and atoms 88 and
131 are the ends of the two hydrophobic chains. The distance from the head to the
middle of the chain ends is 26 Å.

## The numbers

The area per lipid in a fluid POPC bilayer is about 64 Å². With 60 lipids per leaflet,
the membrane has 60 × 64 = 3840 Å², a square with sides of 61.96 Å. Each layer of water
has 15 Å, with 61.96² × 15 = 57600 Å³: at 1 g/mL (29.9 Å³ per molecule), 1925 water
molecules.

The head atoms of the upper leaflet are constrained to be above z = 28 Å, and the chain
ends below z = 6 Å. These planes are 22 Å apart, 85% of the length of the lipid, so the
lipids stay nearly upright (they can tilt by about 30°), with their heads in a narrow
band. Keep these planes close to the length of the lipid: with the heads only 50% of
the length away from the chain ends, the lipids slide along their axes, the chains of
one leaflet push into the other, and the packing gets stuck.

## Command line

```@raw html
<details><summary>The input file <code>bilayer.inp</code> (click to expand)</summary>
```
```@eval
Main.input_file("files/bilayer/bilayer.inp")
```
```@raw html
</details><br>
```

The constraints inside an `atoms ... end atoms` block apply only to the atoms listed:
here, only the head and the chain ends of each lipid are constrained, and the rest of
the molecule follows them. `above plane 0. 0. 1. 28.` is the half-space
`0x + 0y + 1z ≥ 28`, that is, `z ≥ 28`. The chain ends of the two leaflets can meet
anywhere between z = -6 and z = 6, so the leaflets interdigitate a little. The water is
only constrained by one plane, and the box: the water above the membrane fills the box
from z = 30 to its top face, at z = 45, which is in contact, through the periodic
boundary, with the water below the membrane.

### Check the initial placement first

Uncomment the `check_initial_placement` line of the input file, and run it as below
(see [Checking the initial placement](@ref check_initial_placement)): the initial
placement of the molecules is written to the output file (`bilayer.pdb`):

![Initial approximation](images/bilayer_check.webp)

*The initial approximation: the lipids are already oriented by the constraints on their
head and tail atoms, but molecules may overlap.*

Then comment the line again, and run the packing.

### Packing

Run it with:

```
packmol bilayer.inp
```

The packed system (about 28 thousand atoms) is written to `bilayer.pdb`, in about a
minute and a half with 8 threads.

## Julia

The `atom_constraints` keyword of [`structure_type`](@ref) is the equivalent of the
`atoms ... end atoms` blocks: a list of `atom indices => constraints` pairs.

```julia
using Packmol

z = [0., 0., 1.] # the normal of the planes
# Heads (atom 1) above z = 28, chain ends (atoms 88 and 131) below z = 6
upper_leaflet = structure_type("popc.pdb";
    number = 60,
    atom_constraints = [
        [1] => AbovePlane(z, 28.0),
        [88, 131] => BelowPlane(z, 6.0),
    ],
)
# The same, mirrored
lower_leaflet = structure_type("popc.pdb";
    number = 60,
    atom_constraints = [
        [1] => BelowPlane(z, -28.0),
        [88, 131] => AbovePlane(z, -6.0),
    ],
)
water_above = structure_type("water.pdb";
    number = 1925,
    constraints = [AbovePlane(z, 30.0)],
)
water_below = structure_type("water.pdb";
    number = 1925,
    constraints = [BelowPlane(z, -30.0)],
)

system = PackmolSystem([upper_leaflet, lower_leaflet, water_above, water_below];
    output = "bilayer.pdb",
    tolerance = 2.0,
    unitcell = [61.96 0. 0.; 0. 61.96 0.; 0. 0. 90.], # centered at the origin
    # check_initial_placement = true, # uncomment to only place the molecules
)
packmol(system)
```

As with the input file, check the initial placement first, with
`check_initial_placement = true` (commented in the code above).

`BelowPlane(normal, d)` and `AbovePlane(normal, d)` take the parameters of
`below plane` and `above plane`.

!!! tip
    The `Membrane` recipe builds bilayers like this one from a few parameters
    (lipid composition, area per lipid, solvent thickness), with periodic boundary
    conditions. See [Recipes](../recipes/index.md).
