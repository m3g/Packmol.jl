```@meta
CollapsedDocStrings = true
```

# Concentration Unit Conversion

## Overview

These function provide convenient ways to convert between various concentration units using the [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) package for handling physical units. It supports conversions between molarity, molality, mole fraction, mass fraction, volume fraction, and number density.

The primary function for users is `cconvert`. A typical call to `cconvert` takes the concentration value and a `Pair` of strings indicating the input and output units, 
plus the necessary constants (densities, molar masses) necessary for the specific conversion. For example, here we convert a mixture of 40 mL of
ethanol and 60 mL of water (a volume fraction of 40%) to the mole fraction of ethanol:

```@example cconvert0
using Packmol
cconvert(
    40.0,
    "%v/v" => "x";
    M_solute=46.068u"g/mol",   # ethanol
    M_solvent=18.015u"g/mol",  # water
    rho_solute=0.789u"g/mL",   # pure ethanol
    rho_solvent=0.997u"g/mL",  # pure water
)
```

For such conversion the molar masses and the densities of the pure solute and solvent are required. If one of these arguments is not provided,
an error is thrown:

```julia-repl
julia> cconvert(40.0, "%v/v" => "x"; M_solute=46.068u"g/mol", M_solvent=18.015u"g/mol", rho_solute=0.789u"g/mL")
ERROR: UndefKeywordError: keyword argument `rho_solvent` not assigned
```

We can of course convert it back:
```@example cconvert0
cconvert(
    0.17103,
    "x" => "%v/v";
    M_solute=46.068u"g/mol",
    M_solvent=18.015u"g/mol",
    rho_solute=0.789u"g/mL",
    rho_solvent=0.997u"g/mL",
)
```

## Concentration Types

The following types are used internally to dispatch conversion logic but are not typically needed directly by the user when using the string-based `cconvert` method.
To each type a series of aliases is available, which are used as the `String` unit identifiers in the calls to `cconvert`:

* `Molarity`: Represents moles of solute per volume of solution (e.g., mol/L).
* `Molality`: Represents moles of solute per mass of solvent (e.g., mol/kg). 
* `MoleFraction`: Represents moles of solute per total moles (dimensionless fraction, 0-1).
* `MassFraction`: Represents mass of solute per total mass of solution (dimensionless fraction, 0-1). Often expressed as %m/m or w/w. 
* `VolumeFraction`: The volume of the pure solute divided by the sum of the volumes of the pure components, before mixing (dimensionless fraction, 0-1; the [IUPAC definition](https://doi.org/10.1351/goldbook.V06643)). Often expressed as %v/v. Since the volumes of the components are not additive in general, it is not the volume of the solute divided by the volume of the solution: 40 mL of ethanol and 60 mL of water make less than 100 mL of solution, and the volume fraction of ethanol is 0.4.
* `NumberDensity`: Represents number of solute particles per volume of solution (e.g., particles/Å³ or particles/m³). 

The complete list aliases that can be associated to each concentration type is defined in the `Packmol.UNIT_STRINGS`:

```@docs
Packmol.UNIT_STRINGS
```

## API

The main function for performing conversions is `cconvert`.

```@docs
cconvert(::Any, ::Pair{String, String})
```

Convert a concentration `value` from one unit to another using string identifiers. Handles conversion between percentage and fraction representations for dimensionless units.

**Arguments**

* `value`: The numerical value (or `Unitful.Quantity`) of the input concentration.
* `units`: A `Pair` of strings specifying the conversion, e.g., `"mol/L" => "mol/kg"` or `"10 %m/m" => "MoleFraction"`. Supported unit strings include common abbreviations and names (see examples). For dimensionless units (`%m/m`, `%v/v`), providing a "%" sign implies the input `value` is a percentage (e.g., 10.0 for 10%). Otherwise, it's treated as a fraction (0-1). Output format also depends on "%" in the target string.
* `kwargs`: Keyword arguments providing necessary auxiliary data (e.g., `M_solute`, `rho_solution`, `M_solvent`, `rho_solute`). These can be numbers (assuming default units like g/mol, kg/L) or `Unitful.Quantity` objects. The specific required arguments depend on the conversion being performed.

**Returns**

* The converted concentration. This will be a `Unitful.Quantity` for Molarity, Molality, NumberDensity, or a `Real` number for dimensionless fractions (or percentages if requested).

**Examples**

Here, we exemplify conversion of concentration units of Ethanol/Water mixtures.

#### Molarity to Molality (~10 M solution)

```@example cconvert
using Packmol

M_EtOH = 46.068u"g/mol"
M_H2O = 18.015u"g/mol"
rho_EtOH_pure = 0.789u"kg/L" # Density of pure Ethanol
rho_H2O_pure = 0.997u"kg/L" # Density of pure water

rho_sol_10M = 0.90u"kg/L"   # Approx. density for ~10 mol/L EtOH (~58% w/w)
rho_sol_50ww = 0.914u"kg/L" # Approx. density for 50% w/w EtOH
rho_sol_chi02 = 0.94u"kg/L"  # Approx. density for mole fraction χ_EtOH = 0.2

C_in1 = 10.0u"mol/L"
m1 = cconvert(C_in1, "mol/L" => "mol/kg"; M_solute=M_EtOH, rho_solution=rho_sol_10M) # Expect ~22.76 mol/kg
println("10.0 mol/L -> $m1")
```

#### Mass Percent (input %) to Molarity (output Quantity)

```@example cconvert
mass_perc_in = 50.0
C1 = cconvert(mass_perc_in, "%m/m" => "Molarity"; M_solute=M_EtOH, rho_solution=rho_sol_50ww) # Expect ~9.92 mol/L
println("50.0 %m/m -> $C1")
```

#### Mole fraction (input fraction) to Molality (output Quantity)

```@example cconvert
chi_in = 0.2
m2 = cconvert(chi_in, "mole fraction" => "molality"; M_solvent=M_H2O) # Expect ~13.88 mol/kg
println("χ=0.2 -> $m2")
```

#### Molarity (input Quantity) to Volume Percent (output %)
```@example cconvert
C_in2 = 10.0u"M" # Unitful recognizes M as mol/L
vp1 = cconvert(C_in2, "M" => "%v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_10M, rho_solute=rho_EtOH_pure, rho_solvent=rho_H2O_pure)
println("10.0 M -> $vp1 %v/v")
```

#### Number Density (input Quantity) to Molarity
```@example cconvert
N_in = 0.005u"Å^-3"
C2 = cconvert(N_in, "Å^-3" => "Molarity") # Expect ~8.30 mol/L
println("$N_in -> $C2")
```

#### Example using plain numbers for kwargs (assuming default units)
```@example cconvert
m3 = cconvert(0.2, "mole fraction" => "molality"; M_solvent=18.015) # M_H2O in g/mol
println("χ=0.2 (default units) -> $m3")
```

## Required Keyword Arguments for Conversions

The following table indicates the necessary keyword arguments (`kwargs`) required for converting between different concentration units using the `cconvert(value, "UnitIn" => "UnitOut"; kwargs...)` function.

**Keyword Argument Descriptions:**

* `M_solute`: Molar mass of the solute (e.g., in `g/mol`).
* `M_solvent`: Molar mass of the solvent (e.g., in `g/mol`).
* `rho_solution`: Density of the *solution* at the given concentration (e.g., in `g/mL` or `kg/L`).
* `rho_solute`: Density of the *pure solute* (e.g., in `g/mL` or `kg/L`).
* `rho_solvent`: Density of the *pure solvent* (e.g., in `g/mL` or `kg/L`).

| Conversion From => To      | Required Keyword Arguments (`kwargs`)                      | Notes                                           |
| :------------------------- | :--------------------------------------------------------- | :---------------------------------------------- |
| **Molarity =>** |                                                            |                                                 |
| Molarity => Molarity       | *None* | Identity                                        |
| Molarity => Molality       | `M_solute`, `rho_solution`                                 |                                                 |
| Molarity => MoleFraction   | `M_solute`, `M_solvent`, `rho_solution`                    |                                                 |
| Molarity => MassFraction   | `M_solute`, `rho_solution`                                 |                                                 |
| Molarity => VolumeFraction | `M_solute`, `M_solvent`, `rho_solute`, `rho_solvent`, `rho_solution` | Via the mole fraction |
| Molarity => NumberDensity  | *None* | Uses Avogadro constant                          |
| **Molality =>** |                                                            |                                                 |
| Molality => Molarity       | `M_solute`, `rho_solution`                                 |                                                 |
| Molality => Molality       | *None* | Identity                                        |
| Molality => MoleFraction   | `M_solvent`                                                |                                                 |
| Molality => MassFraction   | `M_solute`                                                 |                                                 |
| Molality => VolumeFraction | `M_solute`, `M_solvent`, `rho_solute`, `rho_solvent` | Via the mole fraction |
| Molality => NumberDensity  | `M_solute`, `rho_solution`                                 | Requires intermediate Molarity calc.            |
| **MoleFraction =>** |                                                            |                                                 |
| MoleFraction => Molarity   | `M_solute`, `M_solvent`, `rho_solution`                    |                                                 |
| MoleFraction => Molality   | `M_solvent`                                                | Input `χ` must be < 1                           |
| MoleFraction => MoleFraction | *None* | Identity                                        |
| MoleFraction => MassFraction | `M_solute`, `M_solvent`                                    |                                                 |
| MoleFraction => VolumeFraction | `M_solute`, `M_solvent`, `rho_solute`, `rho_solvent` | |
| MoleFraction => NumberDensity | `M_solute`, `M_solvent`, `rho_solution`                  | Requires intermediate Molarity calc.            |
| **MassFraction =>** |                                                            |                                                 |
| MassFraction => Molarity   | `M_solute`, `rho_solution`                                 |                                                 |
| MassFraction => Molality   | `M_solute`                                                 | Input `mf` must be < 1                          |
| MassFraction => MoleFraction | `M_solute`, `M_solvent`                                    |                                                 |
| MassFraction => MassFraction | *None* | Identity                                        |
| MassFraction => VolumeFraction | `M_solute`, `M_solvent`, `rho_solute`, `rho_solvent` | Via the mole fraction |
| MassFraction => NumberDensity | `M_solute`, `rho_solution`                                 | Requires intermediate Molarity calc.            |
| **VolumeFraction =>** |                                                            |                                                 |
| VolumeFraction => Molarity | `M_solute`, `M_solvent`, `rho_solute`, `rho_solvent`, `rho_solution` | Via the mole fraction |
| VolumeFraction => Molality | `M_solute`, `M_solvent`, `rho_solute`, `rho_solvent` | Via the mole fraction |
| VolumeFraction => MoleFraction | `M_solute`, `M_solvent`, `rho_solute`, `rho_solvent` | |
| VolumeFraction => MassFraction | `M_solute`, `M_solvent`, `rho_solute`, `rho_solvent` | Via the mole fraction |
| VolumeFraction => VolumeFraction | *None* | Identity                                        |
| VolumeFraction => NumberDensity | `M_solute`, `M_solvent`, `rho_solute`, `rho_solvent`, `rho_solution` | Via the mole fraction |
| **NumberDensity =>** |                                                            |                                                 |
| NumberDensity => Molarity   | *None* | Uses Avogadro constant                          |
| NumberDensity => Molality   | `M_solute`, `rho_solution`                                 | Requires intermediate Molarity calc.            |
| NumberDensity => MoleFraction | `M_solute`, `M_solvent`, `rho_solution`                  | Requires intermediate Molarity calc.            |
| NumberDensity => MassFraction | `M_solute`, `rho_solution`                                 | Requires intermediate Molarity calc.            |
| NumberDensity => VolumeFraction | `M_solute`, `M_solvent`, `rho_solute`, `rho_solvent`, `rho_solution` | Via the mole fraction |
| NumberDensity => NumberDensity | *None* | Identity                                        |
