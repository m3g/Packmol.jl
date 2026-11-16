```@meta
CollapsedDocStrings = true
```
# Solute-Solvent system

Here, `SolutionBoxUS` stands for `Solute (U)` and `Solvent (S)`. The density
of the pure solvent is given directly (`g/mL` by default, or `mol/L` for
molarity).

```@setup us_view
using Packmol, PDBTools
test_dir = Packmol.RecipesDirectory * "/test"
system = SolutionBoxUS(
    solute_pdbfile = "$test_dir/data/poly_h.pdb",
    solvent_pdbfile = "$test_dir/data/water.pdb",
    density = 1.0u"g/mL",
)
packed = redirect_stdio(stdout = devnull, stderr = devnull) do
    packmol(system; margin = 20.0, output = joinpath(mktempdir(), "system.pdb"), seed = 1)
end
slab(at) = abs(at.y - packed.unitcell_center[2]) < 5 # a slab of 10 Å through the center
```

```@example us_view
visualize(get_atoms(packed), # hide
    "resname LIG" => (style = :ballandstick,), # hide
    (at -> iswater(at) && slab(at)) => (style = :lines, color = :element, water = true); # hide
    unitcell = UnitCell(packed.unitcell), # hide
    unitcell_origin = packed.unitcell_center - sum(eachcol(packed.unitcell)) / 2, # hide
    rotate = [(:z, 30), (:x, -70)], # hide
) # hide
```

*The example below: the polymer in water, in a cubic periodic box, with only a 10 Å thick section of the water shown, through the center (drag to rotate, scroll to zoom).*

```@docs
SolutionBoxUS
packmol(::SolutionBoxUS)
write_packmol_input(::SolutionBoxUS)
```

## Setting up the system properties

We initialize the system data structure given the PDB files of *one
molecule* of the solute (a polymer, `poly_h.pdb`) and *one molecule* of the
solvent (water):

```@example us
using Packmol
test_dir = Packmol.RecipesDirectory * "/test"
system = SolutionBoxUS(
    solute_pdbfile = "$test_dir/data/poly_h.pdb",
    solvent_pdbfile = "$test_dir/data/water.pdb",
    density = 1.0u"g/mL",
)
```

The molar masses of the solute and solvent can be provided explicitly with
the `solute_molar_mass`/`solvent_molar_mass` keywords. If not, they are
computed from the atom types in the PDB files, which may fail if the mass of
some atom type is unknown.

## Packing the system

The system is built and packed directly, with `packmol`:

```julia
packmol(system; margin = 20.0, output = "system.pdb")
```

This gives the system shown at the top of this page. `output` is the name of the file of
the packed system.

`margin` sets the size of the box from the solute's own bounding box plus
this margin, in every dimension. By default (`pbc = :cubic`) the box is a
cube (all three sides set to the largest of the three margined dimensions);
`pbc = :orthorhombic` instead sizes each side independently, and
`pbc = :dodecahedral`/`pbc = :octahedral` build a rhombic dodecahedron/
truncated octahedron cell, the smallest in which every atom of the solute (centered in the
cell, in the orientation of its file) is at least the margin away from the faces of the
cell (see [The periodic cell](@ref recipes_packing)). For elongated solutes, these cells can be larger
than the cube.

Alternatively, the box size can be given explicitly with
`box_sides = [a, b, c]` (in Å), instead of `margin`.

## Writing a Packmol input file

Instead of packing the system, `write_packmol_input` writes the input file of Packmol
for it (the `input` parameter), to be inspected, edited or run later:

```@example us
write_packmol_input(
    system;
    margin = 20.0,
    input = "box.inp",
    output = "system.pdb",
)
```

The input file is run with the `packmol` executable, or with `packmol("box.inp")` in
Julia (see [Writing a Packmol input file](@ref recipes_input)).
