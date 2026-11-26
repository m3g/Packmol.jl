# Protein in glycerol and water

A protein solvated by a mixture of glycerol and water, 50% in mass, in a periodic box: a
system of the [ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl)
[examples](https://github.com/m3g/ComplexMixturesExamples), where the solvation of the
protein by each component of the mixture is studied. The protein is fixed at the center
of the box, and the mixture fills the rest of it.

```@setup protein_glycerol_view
using Packmol, PDBTools
system = Main.pack_tutorial("protein_glycerol", "protein_glycerol.inp")
```

```@example protein_glycerol_view
visualize(get_atoms(system), # hide
    "protein" => (style = :cartoon,), # hide
    "resname GLYC and y > 0" => (style = :sticks, color = "green"), # hide
    "water and y > 0" => (style = :sticks, color = :element, water = true); # hide
    unitcell = UnitCell(system.unitcell), # hide
    unitcell_origin = system.unitcell_center - sum(eachcol(system.unitcell)) / 2, # hide
    rotate = [(:x, -90)], # hide
) # hide
```

*The protein, glycerol (green) and water, in the half of the box behind the protein, with the edges of the periodic box (drag to rotate, scroll to zoom).*

Files: [`subtilisin.pdb`](https://raw.githubusercontent.com/m3g/MolSimTestFiles/packmol_tutorials-v1/packmol_tutorials/subtilisin.pdb),
[`glycerol.pdb`](https://raw.githubusercontent.com/m3g/MolSimTestFiles/packmol_tutorials-v1/packmol_tutorials/glycerol.pdb),
[`water.pdb`](https://raw.githubusercontent.com/m3g/MolSimTestFiles/packmol_tutorials-v1/packmol_tutorials/water.pdb), and the input file
[`protein_glycerol.inp`](files/protein_glycerol/protein_glycerol.inp).

## The numbers of molecules

The box is a cube with sides of 76 Å, which leaves at least 12 Å of solvent around the
protein (its largest dimension is 51 Å). The mixture of glycerol (92.09 g/mol) and water
(18.02 g/mol), 50% in mass, has a density of 1.126 g/mL. In 1 g of it there are
0.5/92.09 mol of glycerol and 0.5/18.02 mol of water, so its number densities are:

- glycerol: 1.126 × (0.5/92.09) × 0.6022 = 0.00368 molecules per Å³;
- water: 1.126 × (0.5/18.02) × 0.6022 = 0.01882 molecules per Å³.

The protein (27.3 kDa) occupies about 27300/(1.35 × 0.6022) = 33600 Å³ (the density of
proteins is about 1.35 g/mL), so the solution has 76³ − 33600 = 405400 Å³: 1493 molecules of
glycerol and 7630 of water.

!!! tip
    The `SolutionBoxUSC` recipe (see [Recipes](../recipes/solute_solvent_cossolvent_system.md))
    computes these numbers from the density of the mixture as a function of its
    concentration, given as a table, in any concentration unit.

## The protein at the center of the box

```@eval
Main.input_file("files/protein_glycerol/protein_glycerol.inp")
```

The periodic box confines the molecules of the solvent, so they need no constraint. The
protein is fixed with its geometric center (`center`) at the center of the box: all its
atoms must be inside the periodic box (see [Fixed molecules](../constraints/fixed.md)).

### Check the initial placement first

Uncomment the `check_initial_placement` line of the input file, and run it as below
(see [Checking the initial placement](@ref check_initial_placement)): the initial placement,
written to the output file (`protein_glycerol.pdb`) in a fraction of a second, shows the
protein at the center of the box, and the solvent around it. Then comment the line again,
and run the packing.

### Packing

```
packmol protein_glycerol.inp
```

The packing (about 48 thousand atoms) takes about half a second.

## Julia script

```julia
using Packmol

protein = structure_type("subtilisin.pdb";
    number = 1,
    fixed = ([0., 0., 0.], [0., 0., 0.]),
    center = :geometric,
)
glycerol = structure_type("glycerol.pdb"; number = 1493)
water = structure_type("water.pdb"; number = 7630)
system = PackmolSystem([protein, glycerol, water];
    output = "protein_glycerol.pdb",
    tolerance = 2.0,
    unitcell = [76. 0. 0.; 0. 76. 0.; 0. 0. 76.], # centered at the origin
    # check_initial_placement = true, # uncomment to only place the molecules
)
packmol(system)
```

As with the input file, check the initial placement first, with
`check_initial_placement = true` (commented in the code above).
