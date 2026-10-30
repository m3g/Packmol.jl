```@meta
CollapsedDocStrings = true
```
# Coming from MolSimToolkit.jl

The Packmol Input Creator of MolSimToolkit.jl (`MolSimToolkit.PackmolInputCreator`) moved
here. The same systems are built with:

| MolSimToolkit.jl | Packmol.jl |
|---|---|
| `using MolSimToolkit.PackmolInputCreator` | `using Packmol` |
| `PackmolInputCreator.PackmolInputCreatorDirectory` | `Packmol.RecipesDirectory` |
| `SolutionBoxUS(...; density=1.0, density_units="g/mL")` | `SolutionBoxUS(...; density=1.0u"g/mL")` (or `55.5u"mol/L"`; a plain number is in g/mL) |
| `SolutionBoxUSC(...; density_table, concentration_units="x")` | the same |
| `write_packmol_input(system; ..., cubic=true)` | `write_packmol_input(system; ..., pbc=:cubic)` (the default) |
| `write_packmol_input(system; ..., cubic=false)` | `write_packmol_input(system; ..., pbc=:orthorhombic)` |
| `convert_concentration(system, c, "mol/L" => "x")` | the same (molarities are returned with units) |
| `convert_density_table!(system, "mol/L")` | the same; `get_density_table(system)` gives the table in those units |
| `run_packmol("box.inp")` | `packmol("box.inp")`, or `packmol(system; ...)` without an input file |

The concentration units `"mol/L"`, `"x"`, `"vv"` and `"mm"` are accepted, along with
the other aliases of [Concentration Unit Conversion](@ref).
