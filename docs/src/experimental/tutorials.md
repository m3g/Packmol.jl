# Tutorials

!!! warning "Experimental"
    These tutorials use the in-development, pure-Julia packing engine. The
    API shown here can still change.

These tutorials build progressively more complex systems with the
[Julia API](julia_api.md): first a box of water, then a water/urea mixture,
and finally a protein solvated by water and ions. Each packed system is shown
as an interactive 3D view (drag to rotate, scroll to zoom, hover over atoms to
identify them), created with the `visualize` function of
[PDBTools.jl](https://github.com/m3g/PDBTools.jl).

The views are displayed in environments that render HTML: the VS Code plot
pane, Pluto and Jupyter notebooks, and pages like this one. From the REPL,
save a view to a file with `save("view.html", view)` and open it in a web
browser.

## Setup

We load `Packmol` and `PDBTools`, and locate the PDB files of the
molecules used in the tutorials (one molecule per file), which are
distributed with Packmol.jl:

```@example tutorial
using Packmol, PDBTools
data = joinpath(Packmol.RecipesDirectory, "test", "data")
readdir(data)
```

## A box of water

### Water molecules inside a box

The simplest possible system: 1000 water molecules inside a cubic box with
sides of 40 Å. Each kind of molecule is described by a [`structure_type`](@ref),
giving the PDB file of one molecule, the number of copies, and the region where
they must be placed (here, an `InsideBox` defined by its minimum and maximum
corners, centered at the origin):

```@example tutorial
water = structure_type(
    joinpath(data, "water.pdb");
    number = 1000,
    constraints = [InsideBox([-20.0, -20.0, -20.0], [20.0, 20.0, 20.0])],
)
```

The structure types are collected into a [`PackmolSystem`](@ref), which
also holds the global options. The `tolerance` is the minimum distance (in Å)
between atoms of different molecules, and `seed` makes the result
reproducible:

```@example tutorial
system = PackmolSystem([water]; output = "water_box.pdb", tolerance = 2.0, seed = 1)
```

Finally, [`packmol`](@ref) packs the system and writes the output file.
It prints the progress of the optimization (omitted in these tutorials), and
records the outcome in the `status` field of the system, which is
`:packing_ready` if the packing was successful:

```@example tutorial
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
system.status
```

The packed coordinates are retrieved with [`get_atoms`](@ref), which returns
a vector of `PDBTools.Atom`s (the same that would be obtained with
`read_pdb("water_box.pdb")`). We draw the edges of the box with the
`unitcell` keyword of `visualize`, giving the sides of the box and the
position of its corner:

```@example tutorial
atoms = get_atoms(system)
visualize(atoms; unitcell = [40.0, 40.0, 40.0], unitcell_origin = [-20.0, -20.0, -20.0])
```

The molecules fill the box, while no atom of a molecule is closer than 2 Å to
any atom of another. Note that the box is not periodic: molecules at opposite
faces do not see each other, and the density of this box (1000 molecules in
64000 Å³, about 0.47 g/mL) is lower than that of liquid water, leaving room for
the molecules to be easily placed.

### A periodic box at the density of water

For a simulation box we usually want the density of the liquid and periodic
boundary conditions, such that molecules near one face of the box do not
overlap with the periodic images of the molecules at the opposite face.
At 1 g/mL, water has about 0.0334 molecules per Å³, thus 1000 molecules
fit in a cube with sides of:

```@example tutorial
L = round(cbrt(1000 / 0.0334); digits = 2)
```

Periodic boundary conditions are set with the `unitcell` keyword of
`PackmolSystem`, a 3×3 matrix with the cell vectors as columns. As the box
above, the cell is centered at the origin (by default, see `unitcell_center`), and every molecule
is automatically confined to it, so no constraint is needed for the
water molecules:

```@example tutorial
water = structure_type(joinpath(data, "water.pdb"); number = 1000)
system = PackmolSystem(
    [water];
    output = "water_box_pbc.pdb",
    tolerance = 2.0,
    seed = 1,
    unitcell = [L 0 0; 0 L 0; 0 0 L],
)
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
system.status
```

The corner of the cell, used to draw its edges, is at `-L/2` in each direction:

```@example tutorial
atoms = get_atoms(system)
visualize(atoms; unitcell = [L, L, L], unitcell_origin = [-L/2, -L/2, -L/2])
```

The output PDB file contains a `CRYST1` record with the box dimensions, as
required by most simulation packages. See [Periodic boundary
conditions](periodic_boundary_conditions.md) for triclinic, rhombic
dodecahedron, and truncated octahedron cells.

## A mixture of water and urea

Now we build a 4 mol/L solution of urea in water, in a periodic cubic box with
sides of 40 Å. The numbers of molecules follow from the volume of the box,
the concentration of urea, and the density of the solution (about
1.06 g/mL for this concentration). A convenient constant converts the units:
1 mol/L corresponds to 6.022×10⁻⁴ molecules per Å³, and 1 g/mL to 0.6022 g/mol
per Å³:

```@example tutorial
L = 40.0
V = L^3                  # volume of the box, in Å³
c_urea = 4.0             # concentration of urea, in mol/L
ρ = 1.06                 # density of the solution, in g/mL
n_urea = round(Int, c_urea * 6.022e-4 * V)
# the rest of the mass of the box is water
mass_total = ρ * 0.6022 * V                                  # g/mol
mass_urea = n_urea * mass(read_pdb(joinpath(data, "urea.pdb"))) # g/mol
mass_water = mass(read_pdb(joinpath(data, "water.pdb")))
n_water = round(Int, (mass_total - mass_urea) / mass_water)
(; n_urea, n_water)
```

The system has one structure type for each kind of molecule. Both are
confined to the periodic cell:

```@example tutorial
water = structure_type(joinpath(data, "water.pdb"); number = n_water)
urea = structure_type(joinpath(data, "urea.pdb"); number = n_urea)
system = PackmolSystem(
    [water, urea];
    output = "water_urea.pdb",
    tolerance = 2.0,
    seed = 1,
    unitcell = [L 0 0; 0 L 0; 0 0 L],
)
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
system.status
```

By default, `visualize` hides water molecules when other molecules are
present. With `water = true` they are shown as thin lines, while urea is shown
as sticks (with green carbons):

```@example tutorial
atoms = get_atoms(system)
visualize(atoms; water = true, unitcell = [L, L, L], unitcell_origin = [-L/2, -L/2, -L/2])
```

A selection shows only some of the atoms. Here, only the urea molecules, as
balls and sticks:

```@example tutorial
visualize(atoms, "resname URE"; style = :ballandstick, unitcell = [L, L, L], unitcell_origin = [-L/2, -L/2, -L/2])
```

!!! tip
    The [`SolutionBoxUSC`](@ref) recipe computes the numbers of molecules of a
    solute in a mixture of two solvents from the target concentration and a
    table of densities, given in a variety of [concentration
    units](@ref "Concentration Unit Conversion"). See [Recipes](recipes.md).

## A protein solvated by water and ions

The last tutorial builds a typical system for the simulation of a protein: the
protein at the center of a periodic box, surrounded by water, with sodium and
chloride ions at a concentration of 0.15 mol/L, plus the counter-ions required
to neutralize the charge of the protein.

### The protein

We first read the protein structure to obtain its dimensions and its total
charge (computed by PDBTools from the charged residues):

```@example tutorial
protein_file = joinpath(data, "protein.pdb")
protein = read_pdb(protein_file)
protein_charge = round(Int, sum(charge, eachresidue(protein)))
(; natoms = length(protein), maxmin(protein).xlength, protein_charge)
```

The box must be large enough to keep a layer of water around the protein. Here,
we leave at least 10 Å between the protein and each face of the box, such that
the closest distance between the protein and its periodic images is at least
20 Å:

```@example tutorial
margin = 10.0
L = round(maximum(maxmin(protein).xlength) + 2 * margin; digits = 2)
```

### Numbers of molecules

The volume available for the solvent is the volume of the box minus the volume
of the protein, which we estimate assuming that the protein has a density of
1 g/mL. That volume is filled with water at 0.0334 molecules/Å³ and with a salt
concentration of 0.15 mol/L. Since the protein is negatively charged, the
number of sodium ions is increased to neutralize it (and the number of chloride
ions would be increased for a positively charged protein):

```@example tutorial
V_solvent = L^3 - mass(protein) / 0.6022   # Å³
n_water = round(Int, 0.0334 * V_solvent)
n_salt = round(Int, 0.15 * 6.022e-4 * V_solvent)
n_sodium = n_salt + max(0, -protein_charge)
n_chloride = n_salt + max(0, protein_charge)
(; n_water, n_sodium, n_chloride)
```

### Packing

The protein is not moved during the packing: it is a `fixed` structure,
given a position and a rotation (in Euler angles). With `center = :geometric`,
the center of the protein is placed at that position, here the center of the
box (the origin). The water molecules and ions are only confined to the
periodic cell, so they fill the box around the protein:

```@example tutorial
protein_st = structure_type(
    protein_file;
    number = 1,
    fixed = ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
    center = :geometric,
)
water = structure_type(joinpath(data, "water.pdb"); number = n_water)
sodium = structure_type(joinpath(data, "SOD.pdb"); number = n_sodium)
chloride = structure_type(joinpath(data, "CLA.pdb"); number = n_chloride)
system = PackmolSystem(
    [protein_st, water, sodium, chloride];
    output = "solvated_protein.pdb",
    tolerance = 2.0,
    seed = 1,
    unitcell = [L 0 0; 0 L 0; 0 0 L],
)
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
system.status
```

### Visualization

The system has tens of thousands of water molecules, which would hide the
protein. A selection can also be a Julia function of each atom: here we keep
all atoms except the water molecules outside a 6 Å slice through the center of
the box. The protein is shown as a cartoon, colored from the N- (blue) to the
C-terminus (red), and the water molecules as lines:

```@example tutorial
atoms = get_atoms(system)
visualize(atoms, at -> !iswater(at) || abs(at.z) < 3;
    water = true, unitcell = [L, L, L], unitcell_origin = [-L/2, -L/2, -L/2],
)
```

The ions are better seen as spheres. Here we select the protein and the
ions, which are shown in a space-filling representation, with sodium in blue
and chloride in green:

```@example tutorial
visualize(atoms, "protein or resname SOD CLA";
    style = :spheres, unitcell = [L, L, L], unitcell_origin = [-L/2, -L/2, -L/2],
)
```

### The same system with a recipe

The steps above, computing the box size, the charge of the protein, and the
numbers of water molecules and ions, are what the [`SolutionBoxUWI`](@ref)
recipe does. Thus, a system equivalent to the one above is built with:

```@example tutorial
recipe = SolutionBoxUWI(solute_pdbfile = protein_file)
```

```@example tutorial
system = redirect_stdio(stdout = devnull, stderr = devnull) do # hide
system = packmol(recipe; margin = 10.0, ionic_concentration = 0.15u"mol/L", output = "solvated_protein.pdb")
end # hide
system.status
```

The recipe also accepts other ions, other water models, and other shapes of
the periodic cell (rhombic dodecahedron and truncated octahedron), as described
in [Recipes](recipes.md).

## Next steps

- [Constraints](constraints.md): the geometric regions (spheres, cylinders,
  planes, ellipsoids, and more) that molecules can be placed in or kept out
  of, to build interfaces, layers, droplets, and other structures.
- [Periodic boundary conditions](periodic_boundary_conditions.md): triclinic
  and space-saving cell shapes.
- [Recipes](recipes.md): solvation boxes and membranes set up from target
  densities and concentrations.
