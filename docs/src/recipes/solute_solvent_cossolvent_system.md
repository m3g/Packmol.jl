```@meta
CollapsedDocStrings = true
```
# Solute-Solvent-Cossolvent system

Here, `SolutionBoxUSC` stands for `Solute (U)`, `Solvent (S)`, and
`Cossolvent (C)`. Concentrations can be given in molarity (`"mol/L"`), molar
fraction (`"x"`), mass fraction (`"w/w"` or `"mm"`), or volume fraction (`"v/v"` or
`"vv"`); see [Concentration Unit Conversion](@ref) for the full list of unit aliases.

```@setup usc_view
using Packmol, PDBTools
test_dir = Packmol.RecipesDirectory * "/test"
density_table = [
#   x cossolvent (ethanol)     density (g/mL)
             0.0000                 0.9981
             0.0416                 0.9820
             0.0890                 0.9685
             0.1434                 0.9537
             0.2066                 0.9369
             0.2809                 0.9151
             0.3695                 0.8923
             0.4769                 0.8685
             0.6098                 0.8450
             0.7786                 0.8195
             1.0000                 0.7906
]
system = SolutionBoxUSC(
    solute_pdbfile = "$test_dir/data/poly_h.pdb",
    solvent_pdbfile = "$test_dir/data/water.pdb",
    cossolvent_pdbfile = "$test_dir/data/ethanol.pdb",
    density_table = density_table,
    concentration_units = "x",
)
packed = redirect_stdio(stdout = devnull, stderr = devnull) do
    packmol(system; concentration = 0.5, margin = 20.0, output = joinpath(mktempdir(), "system.pdb"), seed = 1)
end
slab(at) = abs(at.y - packed.unitcell_center[2]) < 5 # a slab of 10 Å through the center
```

```@example usc_view
visualize(get_atoms(packed), # hide
    "resname LIG" => (style = :ballandstick,), # hide
    (at -> resname(at) == "ETO" && slab(at)) => (style = :sticks, color = "green"), # hide
    (at -> iswater(at) && slab(at)) => (style = :lines, color = :element, water = true); # hide
    unitcell = UnitCell(packed.unitcell), # hide
    unitcell_origin = packed.unitcell_center - sum(eachcol(packed.unitcell)) / 2, # hide
    rotate = [(:z, 30), (:x, -70)], # hide
) # hide
```

*The example below: the polymer in water and ethanol (green), at a molar fraction of ethanol of 0.5, with only a 10 Å thick section of the solvent shown, through the center (drag to rotate, scroll to zoom).*

```@docs
SolutionBoxUSC
packmol(::SolutionBoxUSC)
write_packmol_input(::SolutionBoxUSC)
convert_concentration
convert_density_table!
get_density_table
```

## Setting up the system properties

The density of the mixture as a function of cossolvent concentration is
given as a table. Here, for a water/ethanol mixture, as a function of the
molar fraction of ethanol:

```@example usc
using Packmol
density_table = [
#   x cossolvent (ethanol)     density (g/mL)
             0.0000                 0.9981
             0.0416                 0.9820
             0.0890                 0.9685
             0.1434                 0.9537
             0.2066                 0.9369
             0.2809                 0.9151
             0.3695                 0.8923
             0.4769                 0.8685
             0.6098                 0.8450
             0.7786                 0.8195
             1.0000                 0.7906
]
```

We then initialize the system data structure, given the PDB files of *one
molecule* of the solute, the solvent (water), and the cossolvent (ethanol):

```@example usc
test_dir = Packmol.RecipesDirectory * "/test"
system = SolutionBoxUSC(
    solute_pdbfile = "$test_dir/data/poly_h.pdb",
    solvent_pdbfile = "$test_dir/data/water.pdb",
    cossolvent_pdbfile = "$test_dir/data/ethanol.pdb",
    density_table = density_table,
    concentration_units = "x", # molar fraction
)
```

As with `SolutionBoxUS`, molar masses are computed from the PDB files unless
given explicitly.

## Packing the system

The system is built and packed directly, at a target cossolvent concentration:

```julia
packmol(system; concentration = 0.5, margin = 20.0, output = "system.pdb")
```

with `concentration` the molar fraction of ethanol (the units of the density table).
This gives the system shown at the top of this page.

`concentration_units` can be passed to interpret `concentration` in units
other than the density table's own (e.g. request a concentration in
`"mol/L"` even though the table above is indexed by molar fraction).
`margin`, `box_sides`, and `pbc` behave as in `SolutionBoxUS`.

## Converting concentrations

`convert_concentration` converts a concentration of the cossolvent between units, using
the density of the solution at that concentration (interpolated from the density table)
and the molar masses of the components. For example, the molarity of ethanol at a molar
fraction of 0.5:

```@example usc
convert_concentration(system, 0.5, "x" => "mol/L")
```

`convert_density_table!` changes the units in which the system takes and shows the
concentrations: after it, `concentration` in `packmol` (and `write_packmol_input`) is read in the new
units, and `get_density_table` returns the table in them:

```@example usc
convert_density_table!(system, "mol/L")
get_density_table(system)
```

!!! note "Volume fractions"
    The volume fraction is the IUPAC one: the volume of the pure cossolvent divided by
    the sum of the volumes of the pure components, before mixing. Computing it requires
    the densities of the pure solvent and cossolvent, so the density table must include
    both (concentrations 0 and 1 in molar fraction).

## Writing a Packmol input file

Instead of packing the system, `write_packmol_input` writes the input file of Packmol
for it, to be inspected, edited or run later (see [Writing a Packmol input file](@ref recipes_input)):

```@example usc
write_packmol_input(
    system;
    concentration = 0.5, # molar fraction of ethanol
    margin = 20.0,
    input = "box.inp",
    output = "system.pdb",
)
```
