# ZIF-8 in water and ethanol

A crystal of the metal–organic framework ZIF-8 in an aqueous solution of ethanol. ZIF-8,
zinc 2-methylimidazolate, has cages of 11.6 Å connected by windows of 3.4 Å, lined by the
methyl groups of the linkers: its pores are hydrophobic, and they take up alcohols but
almost no water, which is why ZIF-8 is studied for the recovery of ethanol and butanol
from water.

The crystal is a slab, periodic in `x` and `y`, and the liquid fills the space between its
two faces, which meet across the periodic boundary in `z`. Ethanol is packed in the pores
and in the liquid; water only in the liquid.

```@setup zif8_final
using Packmol, PDBTools
dir = mktempdir()
for file in readdir(joinpath("files", "zif8"))
    cp(joinpath("files", "zif8", file), joinpath(dir, file))
end
final = redirect_stdio(stdout = devnull, stderr = devnull) do
    packmol(joinpath(dir, "zif8.inp"))
end
```

```@example zif8_final
visualize(get_atoms(final), # hide
    "resname ZIF" => (style = :sticks, color = :element), # hide
    "resname ETOH" => (style = :spheres, color = :element), # hide
    "water" => (style = :sticks, color = :element, water = true); # hide
    unitcell = UnitCell(final.unitcell), # hide
    rotate = [(:z, 30), (:x, -70)], # hide
) # hide
```

*The packed system (drag to rotate, scroll to zoom), built in the [Julia](#Julia) section
below: the framework and water as sticks, ethanol as spheres, colored by element (carbon
gray, hydrogen white, nitrogen blue, oxygen red).*

Files: [`zif8.pdb`](files/zif8/zif8.pdb), [`ethanol.pdb`](files/zif8/ethanol.pdb),
[`water.pdb`](files/zif8/water.pdb), and the input file [`zif8.inp`](files/zif8/zif8.inp).

## The crystal

The structure of ZIF-8 is entry [7111973](https://www.crystallography.net/cod/7111973.html)
of the [Crystallography Open Database](https://www.crystallography.net) (O. Shekhah et al.,
*Chem. Commun.* 50, 2089 (2014),
[doi:10.1039/c3cc47495j](https://doi.org/10.1039/c3cc47495j)): a cubic cell of 17.033 Å,
`Zn₁₂(C₄H₅N₂)₂₄`, with 276 atoms. The full unit cell is obtained from the CIF file with
[Open Babel](https://openbabel.org), as for the [zeolite](zeolite.md):

```
obabel 7111973.cif -O zif8_uc.pdb --fillUC strict
```

The slab, [`zif8.pdb`](files/zif8/zif8.pdb), has 3 × 3 unit cells in `x` and `y` (51.1 Å)
and is 2 unit cells thick in `z`, from `z = 0.5` to `37.7` Å (5364 atoms). It is cut at
whole linkers: the zinc atoms of two unit cells, and every linker bound to them. At each
face, 36 linkers are bound to a single zinc atom, as at the surfaces of real crystals. It was built by the script
`build_slab.jl`, which also keeps one of the two orientations of the methyl groups (the
hydrogen atoms of the methyl groups are disordered in the crystal):

```@raw html
<details><summary>The script <code>build_slab.jl</code> (click to expand)</summary>
```
```@eval
Main.input_file("files/zif8/build_slab.jl")
```
```@raw html
</details><br>
```

The ethanol molecule was built by Open Babel from its SMILES (`obabel -:"CCO" -O
ethanol.pdb --gen3d -h`), and the water molecule is that of the other tutorials.

## The numbers of molecules

The periodic cell is the slab, with 40 Å of liquid above it: `pbc 51.099 51.099 80.0`
(from the origin to these sides).

**The liquid** is the region 2 Å beyond the faces of the crystal (`39.66 < z < 80`),
about 105000 Å³. An aqueous solution of 20% ethanol in volume has 217 molecules of
ethanol (0.789 g/mL) and 2808 of water (0.997 g/mL) in it.

**The pores**: the volume of the crystal available to ethanol is given by
[`free_volume`](@ref), which samples the region of each structure type and discounts the
space taken by the fixed molecules (see [Volume estimates](@ref)):

```
packmol --free-volume zif8.inp
```
```
    Structure                 Number     Volume (Å³)     ± error  Density (g/mL)
    ethanol.pdb                  405           39308        48.7           0.788
    ethanol.pdb                  217          104861         7.0           0.158
    water.pdb                   2808          104854         7.1           0.801
```

The crystal has 39300 Å³ available to ethanol, 0.48 cm³ per gram of ZIF-8. Filled with
liquid ethanol, at 0.789 g/mL, it holds 405 molecules, 22 per unit cell. This is an upper
bound of what the crystal adsorbs: the adsorbed phase is less dense than the liquid. The
other two lines check the liquid: ethanol and water fill it at 0.158 + 0.801 = 0.959
g/mL, the density of the mixture.

## Command line

```@raw html
<details><summary>The input file <code>zif8.inp</code> (click to expand)</summary>
```
```@eval
Main.input_file("files/zif8/zif8.inp")
```
```@raw html
</details><br>
```

The crystal is a fixed molecule (see [Fixed molecules](../constraints/fixed.md)), with its
atoms inside the periodic cell. Ethanol has two structure blocks: one in the region of the
crystal, where the distances to the atoms of the framework put the molecules in the
cages, and one in the liquid. Water has a single block, in the liquid: it is the
constraint that keeps water out of the pores. The `radius` of the atoms of the crystal,
1.5 Å, keeps the atoms of the guests 2.5 Å from those of the framework.

### Check the initial placement first

Uncomment the `check_initial_placement` line of the input file, and run it as below
(see [Checking the initial placement](@ref check_initial_placement)): the initial
placement, written to the output file in a fraction of a second, already has the ethanol
molecules in the cages, and the liquid above the crystal.

When it does (see the views in the [Julia](#Julia) section), comment the line again, and
run the packing.

### Packing

```
packmol zif8.inp
```

The packing (19386 atoms) takes half a second.

## Julia

```@setup zif8
for file in readdir(joinpath("files", "zif8"))
    cp(joinpath("files", "zif8", file), file; force = true)
end
```

The crystal, fixed, and the ethanol in its pores:

```@example zif8
using Packmol, PDBTools

zif8 = structure_type("zif8.pdb";
    number = 1,
    fixed = ([0., 0., 0.], [0., 0., 0.]),
    radius = 1.5,
)
pores = structure_type("ethanol.pdb";
    number = 405,
    constraints = [InsideBox([0., 0., 0.], [51.099, 51.099, 38.16])],
)
nothing # hide
```

The liquid, a region with ethanol and water:

```@example zif8
liquid = InsideBox([0., 0., 39.66], [51.099, 51.099, 80.0])
ethanol = structure_type("ethanol.pdb"; number = 217, constraints = [liquid])
water = structure_type("water.pdb"; number = 2808, constraints = [liquid])
nothing # hide
```

The system, in the periodic cell of `pbc 51.099 51.099 80.0` (from the origin), first
with only the initial placement:

```@example zif8
system = PackmolSystem([zif8, pores, ethanol, water];
    output = "zif8_ethanol_water.pdb",
    tolerance = 2.0,
    unitcell = [51.099 0. 0.; 0. 51.099 0.; 0. 0. 80.0],
    unitcell_center = [25.5495, 25.5495, 40.0],
    check_initial_placement = true, # first, only the initial placement
)
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
system.status
```

The initial placement, with the `visualize` function of
[PDBTools.jl](https://github.com/m3g/PDBTools.jl) (drag to rotate, scroll to zoom): the
ethanol molecules are already in the cages, and the liquid above the crystal:

```@example zif8
function show_system(system)
    return visualize(get_atoms(system),
        "resname ZIF" => (style = :sticks, color = :element),
        "resname ETOH" => (style = :spheres, color = :element),
        "water" => (style = :sticks, color = :element, water = true);
        unitcell = UnitCell(system.unitcell),
        rotate = [(:z, 30), (:x, -70)],
    )
end
show_system(system)
```

Then the packing:

```@example zif8
system.check_initial_placement = false
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
show_system(system)
```

`free_volume(system)` gives the volumes of the table above.
