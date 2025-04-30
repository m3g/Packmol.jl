using Unitful
export @u_str 
@unit molecule "molecule" molecule 1u"mol" / 6.02e23 false

Unitful.register(@__MODULE__)
function __init__()
    return Unitful.register(@__MODULE__)
end  

export cconvert

const unit_names = Dict{String,String}(
    "mol/L" => "molarity",
    "mol/kg" => "molality",
    "x" => "molar fraction",
    "vv" => "volume fraction",
    "mm" => "mass fraction",
    "molecule/Å^3" => "molecule/Å^3",
    "molecule/angstrom^3" => "molecule/Å^3"
)
unit_name(u::String) = unit_names[u]

abstract type AbstractUnit end
struct Molarity <: AbstractUnit value::Float64 end
struct Molality <: AbstractUnit value::Float64 end
struct MolarFraction <: AbstractUnit value::Float64 end
struct VolumeFraction <: AbstractUnit value::Float64 end
struct MassFraction <: AbstractUnit value:: Float64 end
struct MoleculeAngs3 <: AbstractUnit value::Float64 end

const unit_type = Dict{String,DataType}(
    "mol/L" => Molarity,
    "mol/kg" => Molality,
    "x" => MolarFraction,
    "vv" => VolumeFraction,
    "mm" => MassFraction,
    # molecular-length units
    "molecule/A^3" => MoleculeAngs3,
    "molecule/Å^3" => MoleculeAngs3,
    "molecule/angs^3" => MoleculeAngs3,
    "molecule/angstrom^3" => MoleculeAngs3,
)

"""
    cconvert(value::Number, units::Pair{String,String}; kargs...)

Convert a number `value` from one unit to another. The units are specified as a pair of strings, where the first string 
is the source unit and the second string is the target unit. 

Each specific conversion requires additional keyword arguments, which are passed as `kargs...`.
The function will throw an error if the conversion is not possible or if the required keyword arguments are not provided.

`Unitful` is units are supported, which can help with the proper conversion of units.

# Example:

```jldoctest
julia> using Packmol

julia> cconvert(55.56u"mol/L", "mol/L" => "mol/kg"; density=1.0u"g/mL") # pure water
55.56 mol kg^-1

# Convert from molarity to volume fraction (ethanol-water mixture)
julia> cconvert(0.5u"mol/L", "mol/L" => "vv"; 
                  density=0.993u"g/mL",
                  density_pure_solute=0.789u"g/mL", # pure ethanol
                  density_pure_solvent=1.0u"g/mL", # pure water
                  solute_molar_mass=46.07u"g/mol", # ethanol
              )
0.029219722974320184

# Convert back from volume fraction to molarity
julia> cconvert(ans, "vv" => "mol/L"; 
           density=0.993u"g/mL", 
           density_pure_solute=0.789u"g/mL", 
           density_pure_solvent=1.0u"g/mL", 
           solute_molar_mass=46.07u"g/mol"
        )
0.5 mol L^-1
```

"""
cconvert(x::Number, units::Pair{String,String}; kargs...) = cconvert(x, unit_type[units[1]], unit_type[units[2]]; kargs...)

#
# Convert Molarity to other units
#
cconvert(x::Number, ::Type{Molarity}, ::Type{Molarity}) = x
cconvert(x::Number, ::Type{Molarity}, ::Type{Molality}; density) = upreferred(x/density) # mol/L * mL/g = mol/kg 
function cconvert(x::Number, ::Type{Molarity}, ::Type{MolarFraction}; 
    density_pure_solvent, solute_molar_mass, solvent_molar_mass
)
    ns = (density_pure_solvent - x * solute_molar_mass) / solvent_molar_mass
    return upreferred(x / (x + ns))
end
function cconvert(x::Number, ::Type{Molarity}, ::Type{VolumeFraction}; 
    density, density_pure_solute, density_pure_solvent, solute_molar_mass
)
    mu = x * solute_molar_mass
    vu = mu / density_pure_solute
    ms = density - mu 
    vs = ms / density_pure_solvent
    return upreferred(vu / (vu + vs))
end
function cconvert(x::Number, ::Type{Molarity}, ::Type{MassFraction}; 
    density, solute_molar_mass
)
    mu = x * solute_molar_mass
    ms = density - mu
    return upreferred(mu / (mu + ms))
end

#
# Convert from other units to Molarity
#
cconvert(x::Number, ::Type{Molality}, ::Type{Molarity}; density) = 
    x isa Quantity ? uconvert(u"mol/L", x * density) : x * density # mol/kg * g/mL = mol/L
function cconvert(x::Number, ::Type{MolarFraction}, ::Type{Molarity}; 
    density, solute_molar_mass, solvent_molar_mass
)
    mol_per_g = inv(x * solute_molar_mass + (1 - x) * solvent_molar_mass)
    mol_per_ml = mol_per_g * density
    mol_of_solute_per_ml = x * mol_per_ml
    return density isa Quantity ? uconvert(u"mol/L", mol_of_solute_per_ml) : mol_of_solute_per_ml * 1e3
end
function cconvert(x::Number, ::Type{VolumeFraction}, ::Type{Molarity}; 
    density, density_pure_solute, density_pure_solvent, solute_molar_mass
)
    v_u = x
    v_s = 1 - x
    m_u = v_u * density_pure_solute # in 1 mL
    n_u = m_u / solute_molar_mass
    m_s = v_s * density_pure_solvent # in 1 mL
    tot_mass = m_u + m_s
    tot_vol = tot_mass / density
    return density isa Quantity ? uconvert(u"mol/L", n_u / tot_vol) : n_u / tot_vol * 1e3
end

#
# Convert from molarity to molecules/Å^3
#
cconvert(x::Number, ::Type{Molarity}, ::Type{MoleculeAngs3}) = uconvert(u"molecule/Å^3", x)

@testitem "cconvert" begin
    using Packmol, Unitful
    # Molarity to all
    @test cconvert(1.0u"mol/L", "mol/L" => "mol/L") == 1.0u"mol/L"
    @test cconvert(1.0u"mol/L", "mol/L" => "mol/kg"; density=1.0u"g/mL") == 1.0u"mol/kg"
    @test cconvert(1.0u"mol/L", "mol/L" => "mol/kg"; density=0.5u"g/ml") == 2.0u"mol/kg"
    @test_throws UndefKeywordError cconvert(1.0, "mol/L" => "mol/kg")
    @test cconvert(55.56u"mol/L", "mol/L" => "x"; density_pure_solvent=1.0u"g/mL", solute_molar_mass=18u"g/mol", solvent_molar_mass=18u"g/mol") ≈ 1 atol=1e-3
    @test cconvert(55.56u"mol/L", "mol/L" => "vv"; density_pure_solvent=1.0u"g/mL", density_pure_solute=1.0u"g/mL", solute_molar_mass=18u"g/mol", density=1.0u"g/mL") ≈ 1 atol=1e-3
    @test cconvert(55.56u"mol/L", "mol/L" => "mm"; density=1.0u"g/mL", solute_molar_mass=18u"g/mol") ≈ 1 atol=1e-3
    @test cconvert(55.56u"mol/L", "mol/L" => "molecule/Å^3") ≈ 0.0334u"molecule/Å^3" atol=1e-3u"molecule/Å^3"
    # All to molarity
    @test cconvert(1u"mol/kg", "mol/kg" => "mol/L"; density=1.00u"g/mL") == 1u"mol/L"
    @test cconvert(1u"mol/kg", "mol/kg" => "mol/L"; density=2.00u"g/mL") == 2.0u"mol/L"
    @test_throws UndefKeywordError cconvert(1.0, "mol/kg" => "mol/L")
    @test cconvert(1, "x" => "mol/L"; density=1u"g/ml", solute_molar_mass=18u"g/mol", solvent_molar_mass=18u"g/mol") ≈ 55.56u"mol/L" atol=1e-2u"mol/L"
    @test cconvert(1, "x" => "mol/L"; density=1, solute_molar_mass=18, solvent_molar_mass=18) ≈ 55.56 atol=1e-2
    @test_throws UndefKeywordError cconvert(1, "x" => "mol/L"; solute_molar_mass=18u"g/mol", solvent_molar_mass=18u"g/mol")
    @test cconvert(1, "vv" => "mol/L"; density=1u"g/mL", density_pure_solute=1u"g/mL", density_pure_solvent=1u"g/mL", solute_molar_mass=18u"g/mol") ≈ 55.56u"mol/L" atol=1e-2u"mol/L"
    @test cconvert(1, "vv" => "mol/L"; density=1, density_pure_solute=1, density_pure_solvent=1, solute_molar_mass=18) ≈ 55.56 atol=1e-2

end

#=
    interpolate_concentration(system, x)

Obtain by interpolation the value of of ρ[:,2] that corresponds
to a the best estimate given a value of x corresponding to 
to the domain ρ[:,1].

=#
function interpolate_concentration(system, x)
    ρ = system.density_table
    i = findfirst(d -> d >= x, @view(ρ[:, 1]))
    i == firstindex(@view(ρ[:, 1])) && return ρ[i, 2]
    ρ[i] == x && return ρ[i, 2]
    dρdx = (ρ[i, 2] - ρ[i-1, 2]) / (ρ[i, 1] - ρ[i-1, 1])
    d = ρ[i-1, 2] + dρdx * (x - ρ[i-1, 1])
    return d
end

# Adjust x to be in the [0,1] range, to avoid problems with numerical precision
fixrange(x) = x < 0 ? 0 : (x > 1 ? 1 : x)

# System types
#include("./SolutionBoxUS.jl")
#include("./SolutionBoxUSC.jl")
#include("./SolutionBoxUWI.jl")

# Tests
#include("./test/runtests.jl")
