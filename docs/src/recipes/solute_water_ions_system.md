```@meta
CollapsedDocStrings = true
```
# Solute-Water-Ions system

Here, `SolutionBoxUWI` stands for `Solute (U)`, `Water (W)`, and `Ions (I)`:
a solute solvated by water and a background salt at a target ionic
concentration (`0.15u"mol/L"`, physiological saline, by default). The
solute's own net charge (`solute_charge`, computed from the PDB file's
`charge` records unless given explicitly) is automatically neutralized with
extra counter-ions, on top of the bulk salt.

```@setup uwi_view
using Packmol, PDBTools
test_dir = Packmol.RecipesDirectory * "/test"
system = SolutionBoxUWI(
    solute_pdbfile = "$test_dir/data/poly_h.pdb",
    solute_charge = 2,
)
packed = redirect_stdio(stdout = devnull, stderr = devnull) do
    packmol(system; ionic_concentration = 0.15u"mol/L", margin = 20.0, output = joinpath(mktempdir(), "system.pdb"), seed = 1)
end
```

```@example uwi_view
visualize(get_atoms(packed), # hide
    "resname LIG" => (style = :ballandstick,), # hide
    "resname SOD" => (style = :spheres, color = "orange"), # hide
    "resname CLA" => (style = :spheres, color = "green"), # hide
    (at -> iswater(at) && element(at) == "O" && resnum(at) % 10 == 0) => (style = :spheres, color = "lightblue", opacity = 0.6, water = true); # hide
    unitcell = UnitCell(packed.unitcell), # hide
    unitcell_origin = packed.unitcell_center - sum(eachcol(packed.unitcell)) / 2, # hide
    rotate = [(:z, 30), (:x, -70)], # hide
) # hide
```

*The example below: the polymer in water, with sodium (orange) and chloride (green) ions, and the oxygen atoms of one in ten water molecules (light blue) (drag to rotate, scroll to zoom).*

```@docs
SolutionBoxUWI
packmol(::SolutionBoxUWI)
write_packmol_input(::SolutionBoxUWI)
```

## Setting up the system properties

By default the ions are sodium and chloride (charges `+1`/`-1`), generated
automatically, and the solution density as a function of ionic strength is
looked up from a table of aqueous NaCl densities:

```@example uwi
using Packmol
test_dir = Packmol.RecipesDirectory * "/test"
system = SolutionBoxUWI(
    solute_pdbfile = "$test_dir/data/poly_h.pdb",
    solute_charge = 2,
)
```

Custom cation/anion PDB files can be given with the `cation_pdbfile`/
`anion_pdbfile` keywords (with their `cation_charge`/`anion_charge`); a
matching `density_table` should then be provided too, since the default one
is specific to NaCl.

## Packing the system

The system is built and packed directly, at a target ionic concentration:

```julia
packmol(system;
    ionic_concentration = 0.15u"mol/L",
    margin = 20.0,
    output = "system.pdb",
)
```

This gives the system shown at the top of this page.

Enough cations or anions (of a single sign) are added on top of the bulk
salt to exactly neutralize `system.solute_charge`; a `SolutionBoxUWI` with
default (monovalent) ions can neutralize any integer solute charge. An
`ArgumentError` is raised if the requested charge cannot be exactly
neutralized with the ions' charges (for instance, a solute charge of `-3`
with only doubly-charged cations available). `margin`, `box_sides`, and
`pbc` behave as in `SolutionBoxUS`.

## Writing a Packmol input file

Instead of packing the system, `write_packmol_input` writes the input file of Packmol
for it, to be inspected, edited or run later (see [Writing a Packmol input file](@ref recipes_input)):

```@example uwi
write_packmol_input(
    system;
    ionic_concentration = 0.15u"mol/L",
    margin = 20.0,
    input = "box.inp",
    output = "system.pdb",
)
```
