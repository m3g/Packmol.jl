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

```@docs
SolutionBoxUWI
write_packmol_input(::SolutionBoxUWI)
packmol(::SolutionBoxUWI)
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

Finally, we generate the input file at a target ionic concentration:

```@example uwi
write_packmol_input(
    system;
    ionic_concentration = 0.15u"mol/L",
    margin = 20.0,
    input = "box.inp",
    output = "system.pdb",
)
```

or, equivalently, pack it directly with `packmol(system; ionic_concentration=0.15u"mol/L", margin=20.0, output="system.pdb")`.

Enough cations or anions (of a single sign) are added on top of the bulk
salt to exactly neutralize `system.solute_charge`; a `SolutionBoxUWI` with
default (monovalent) ions can neutralize any integer solute charge. An
`ArgumentError` is raised if the requested charge cannot be exactly
neutralized with the ions' charges (for instance, a solute charge of `-3`
with only doubly-charged cations available). `margin`, `box_sides`, and
`pbc` behave as in `SolutionBoxUS`.
