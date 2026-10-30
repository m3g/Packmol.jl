```@meta
CollapsedDocStrings = true
```
# Solute-Solvent system

Here, `SolutionBoxUS` stands for `Solute (U)` and `Solvent (S)`. The density
of the pure solvent is given directly (`g/mL` by default, or `mol/L` for
molarity).

```@docs
SolutionBoxUS
write_packmol_input(::SolutionBoxUS)
packmol(::SolutionBoxUS)
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

Finally, we generate an input file for Packmol with:

```@example us
write_packmol_input(
    system;
    margin = 20.0,
    input = "box.inp",
    output = "system.pdb",
)
```

or build and pack the system directly, with `packmol`:

```julia
packmol(system; margin = 20.0, output = "system.pdb")
```

The `input` parameter (`write_packmol_input` only) is the name of the
generated Packmol input file, and `output` is the name assigned to the
packed system.

`margin` sets the size of the box from the solute's own bounding box plus
this margin, in every dimension. By default (`pbc = :cubic`) the box is a
cube (all three sides set to the largest of the three margined dimensions);
`pbc = :orthorhombic` instead sizes each side independently, and
`pbc = :dodecahedral`/`pbc = :octahedral` build a rhombic dodecahedron/
truncated octahedron cell of that same largest side (see [Running
Packmol](@ref) above).

Alternatively, the box size can be given explicitly with
`box_sides = [a, b, c]` (in Å), instead of `margin`.
