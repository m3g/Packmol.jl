```@meta
CollapsedDocStrings = true
```
# Membrane system

`Membrane` builds a lipid bilayer or monolayer under orthorhombic periodic
boundary conditions, with the membrane normal along z. Unlike the solvation
recipes above, every parameter is given up front in the `Membrane`
constructor itself — `write_packmol_input`/`packmol` take only `input`/
`output` (and, for `packmol`, any packing-engine keyword).

```@docs
Membrane
write_packmol_input(::Membrane)
packmol(::Membrane)
```

## Setting up the system properties

Each lipid is identified by a single "head" atom and a single "tail" atom
(1-based indices into its own PDB file); the Euclidean distance between them
sets that lipid's head-to-tail length, and the *longest* one (across all
lipids given) sets the membrane's nominal thickness. Here, a single lipid
type, solvated by water:

```@example membrane
using Packmol
test_dir = Packmol.RecipesDirectory * "/test"
system = Membrane(
    type = :bilayer,
    lipids = "$test_dir/data/lipid.pdb",
    lipid_head = [31],  # the head-group ring
    lipid_tail = [1],   # the far end of the acyl chain
    lipid_molar_ratio = [1.0],
    area_per_lipid = 60.0, # Å², xor total_area
    total_lipids = 200,    # xor total_area, counted over both leaflets
    solvent = "$test_dir/data/water.pdb",
    solvent_layer_width = 20.0, # Å
    solvent_density = 1.0u"g/mL",
)
```

Mixing lipids or solvents just means passing `Vector`s of the same length to
`lipids`/`lipid_head`/`lipid_tail`/`lipid_molar_ratio` (and, respectively,
`solvent`/`solvent_molar_ratio`) instead of scalars/single files — see the
`Membrane` docstring above for the full parameter list, including
`flexibility` (how much tilt/wobble each lipid is allowed while its head and
tail are pinned near their target planes) and the `total_area`/
`total_lipids` either-or.

Finally, generate the input file, or pack it directly:

```@example membrane
write_packmol_input(system; input = "membrane.inp", output = "membrane.pdb")
```

```julia
packmol(system; output = "membrane.pdb")
```
