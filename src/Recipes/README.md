# Recipes

Systems set up from densities and concentrations instead of numbers of molecules. Each
recipe is a type (a subtype of `Recipe`, built with keyword arguments only) with a
`write_packmol_input` method (writes an input file) and `packmol(recipe)` (builds and
packs the system, with no file left behind). The recipes are part of the `Packmol` module.

| File | Contents |
|---|---|
| `Recipes.jl` | `Recipe`, the checks of the keyword arguments, `packmol(recipe)`, and the includes. |
| `SolutionBoxUS.jl` | A solute in a solvent. |
| `SolutionBoxUSC.jl` | A solute in a mixture of a solvent and a cossolvent, at a concentration of the cossolvent; the conversion of its concentrations and of its density table between units (`convert_concentration`, `convert_density_table!`, `get_density_table`). |
| `SolutionBoxUWI.jl` | A solute in water with ions, at a salt concentration, neutralizing the charge of the solute. |
| `Membrane.jl` | A lipid bilayer or monolayer, of one or more lipids, solvated, in an orthorhombic cell. |
| `MembraneProtein.jl` | A membrane protein (oriented as in OPM) in a bilayer, with water and ions. |
| `DensityTable.jl` | Tables of densities as a function of the concentration. |
| `concentration_units.jl` | Concentration units (molarity, molality, mole and mass fractions, ...) and the conversions between them, given the densities of the solutions (`cconvert`). |
| `lipids/` | The structures of the built-in lipids (see its README). |
| `test/` | Data files used by the examples of the documentation (not test code). |
