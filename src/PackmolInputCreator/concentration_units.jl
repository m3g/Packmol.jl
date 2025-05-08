using Unitful
using OrderedCollections: OrderedDict

export @u_str
export cconvert # Export the user-friendly string-based function


# ==============================================================================
# Define Concentration Unit Types (as Dispatch Tags)
# ==============================================================================
abstract type AbstractConcentrationUnit end

struct Molarity <: AbstractConcentrationUnit end        # Represents mol/Volume (e.g., mol/L)
struct Molality <: AbstractConcentrationUnit end        # Represents mol/Mass (e.g., mol/kg solvent)
struct MoleFraction <: AbstractConcentrationUnit end    # Represents mol/mol (dimensionless fraction 0-1)
struct MassFraction <: AbstractConcentrationUnit end    # Represents mass/mass (dimensionless fraction 0-1)
struct VolumeFraction <: AbstractConcentrationUnit end  # Represents volume/volume (dimensionless fraction 0-1)
struct NumberDensity <: AbstractConcentrationUnit end   # Represents count/Volume (e.g., Å^-3 or m^-3)

# ==============================================================================
# Default Units and Helper Function
# ==============================================================================

_oneunit(::Type{Molarity}) = 1u"mol/L"
_oneunit(::Type{Molality}) = 1u"mol/kg"
_oneunit(::Type{MoleFraction}) = 1
_oneunit(::Type{MassFraction}) = 1
_oneunit(::Type{VolumeFraction}) = 1
_oneunit(::Type{NumberDensity}) = 1u"Å^-3"

const DEFAULT_MOLARITY_UNIT = u"mol/L"
const DEFAULT_MOLALITY_UNIT = u"mol/kg"
const DEFAULT_MASS_UNIT = u"g/mol" # For molar mass
const DEFAULT_DENSITY_UNIT = u"g/mL" 
const DEFAULT_VOLUME_UNIT = u"L"
const DEFAULT_MASS_QUANTITY_UNIT = u"kg" # For mass of solvent/solution
const DEFAULT_NUMBERDENSITY_UNIT = u"Å^-3"
const AVOGADRO = Unitful.Na # Avogadro constant

#=
    _ensure_unit(value, default_unit)

Helper to attach a default unit if the input is a plain number,
or convert to the default unit if it's already a Quantity.
=#
function _ensure_unit(val, default_unit::Unitful.Units)
    if val isa Quantity
        return uconvert(default_unit, val)
    elseif val isa Real
        return val * default_unit
    else
        throw(ArgumentError("Input must be a Real number or a Unitful.Quantity, got $(typeof(val))"))
    end
end

#=
    _ensure_unit(value, default_unit::Unitful.Quantity)

Method to handle cases where the default is already a quantity (like Avogadro).
=#
function _ensure_unit(val, default_quantity::Unitful.Quantity)
     if val isa Quantity
         # This case is ambiguous - should we compare dimensions or just return val?
         # Let's assume the user provided the quantity correctly.
         return val
     elseif val isa Real
         # This is also ambiguous - multiply the number by the default quantity?
         # Unlikely to be intended. Let's error.
         throw(ArgumentError("Cannot apply default quantity $default_quantity to plain number $val."))
     else
        throw(ArgumentError("Input must be a Real number or a Unitful.Quantity, got $(typeof(val))"))
    end
end


# ==============================================================================
# Core Conversion Logic using Multiple Dispatch
# ==============================================================================

# --- Identity Conversions ---
function cconvert(x, ::Type{U}, ::Type{U}; kwargs...) where {U<:AbstractConcentrationUnit}
    # Ensure consistent return type (Quantity or Real for fractions)
     if U == Molarity return _ensure_unit(x, DEFAULT_MOLARITY_UNIT) end
     if U == Molality return _ensure_unit(x, DEFAULT_MOLALITY_UNIT) end
     if U == NumberDensity return _ensure_unit(x, DEFAULT_NUMBERDENSITY_UNIT) end
     if U in (MoleFraction, MassFraction, VolumeFraction)
        # Expect dimensionless fraction (Real or DimensionlessQuantity)
        return x isa Quantity ? ustrip(NoUnits, x) : x
     end
     # Fallback - should not be reached if all types are handled
     return x
end

# --- Molarity Conversions ---
function cconvert(
    C_val, ::Type{Molarity}, ::Type{Molality};
    M_solute, rho_solution, kwargs...
)
    C = _ensure_unit(C_val, DEFAULT_MOLARITY_UNIT)
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)
    rho_s = _ensure_unit(rho_solution, DEFAULT_DENSITY_UNIT)

    V0 = 1u"L" # Basis volume
    n_solute = uconvert(u"mol", C * V0)
    m_solution = uconvert(DEFAULT_MASS_QUANTITY_UNIT, rho_s * V0)
    m_solute = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solute * Ms)

    if m_solute >= m_solution throw(ArgumentError("Solute mass ($m_solute) >= solution mass ($m_solution)")) end
    m_solvent = m_solution - m_solute
    if m_solvent <= 0 * unit(m_solvent) throw(ArgumentError("Non-positive solvent mass ($m_solvent)")) end

    molality = n_solute / m_solvent
    return uconvert(DEFAULT_MOLALITY_UNIT, molality)
end

function cconvert(
    C_val, ::Type{Molarity}, ::Type{MoleFraction};
    M_solute, M_solvent, rho_solution, kwargs...
)
    C = _ensure_unit(C_val, DEFAULT_MOLARITY_UNIT)
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)
    Msv = _ensure_unit(M_solvent, DEFAULT_MASS_UNIT)
    rho_s = _ensure_unit(rho_solution, DEFAULT_DENSITY_UNIT)

    V0 = 1u"L" # Basis volume
    n_solute = uconvert(u"mol", C * V0)
    m_solution = uconvert(DEFAULT_MASS_QUANTITY_UNIT, rho_s * V0)
    m_solute = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solute * Ms)

    if m_solute >= m_solution throw(ArgumentError("Solute mass ($m_solute) >= solution mass ($m_solution)")) end
    m_solvent = m_solution - m_solute
    if m_solvent <= 0 * unit(m_solvent) throw(ArgumentError("Non-positive solvent mass ($m_solvent)")) end

    n_solvent = uconvert(u"mol", m_solvent / Msv)
    n_total = n_solute + n_solvent

    return n_total == 0 * unit(n_total) ? 0.0 : ustrip(NoUnits, n_solute / n_total)
end

function cconvert(
    C_val, ::Type{Molarity}, ::Type{MassFraction};
    M_solute, rho_solution, kwargs...
)
    C = _ensure_unit(C_val, DEFAULT_MOLARITY_UNIT)
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)
    rho_s = _ensure_unit(rho_solution, DEFAULT_DENSITY_UNIT)

    V0 = 1u"L" # Basis volume
    n_solute = uconvert(u"mol", C * V0)
    m_solution = uconvert(DEFAULT_MASS_QUANTITY_UNIT, rho_s * V0)
    m_solute = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solute * Ms)

    if m_solution <= 0 * unit(m_solution) throw(ArgumentError("Non-positive solution mass ($m_solution)")) end
    # Allow m_solute == m_solution for pure solute case if C is very high? No, rho_s would be rho_solute.
    if m_solute > m_solution throw(ArgumentError("Solute mass ($m_solute) > solution mass ($m_solution)")) end

    mf = m_solute / m_solution # Should be dimensionless
    return ustrip(NoUnits, mf)
end

function cconvert(
    C_val, ::Type{Molarity}, ::Type{VolumeFraction};
    M_solute, rho_solute, kwargs...
)
    # This conversion assumes the definition V_solute / V_solution
    # V_solution is the basis (e.g. 1L). V_solute depends on density of PURE solute.
    # rho_solution is NOT directly needed here if we assume C refers to the final volume V0.
    C = _ensure_unit(C_val, DEFAULT_MOLARITY_UNIT)
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)
    rho_pure_s = _ensure_unit(rho_solute, DEFAULT_DENSITY_UNIT)

    V0 = 1u"L" # Basis volume (final solution volume)
    n_solute = uconvert(u"mol", C * V0)
    m_solute = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solute * Ms)
    V_solute = uconvert(DEFAULT_VOLUME_UNIT, m_solute / rho_pure_s)

    # Volume fraction assumes additivity, or more commonly relates V_solute to V_solution
    vf = V_solute / V0 # Should be dimensionless
    return ustrip(NoUnits, vf)
end

function cconvert(C_val, ::Type{Molarity}, ::Type{NumberDensity}; kwargs...)
    C = _ensure_unit(C_val, DEFAULT_MOLARITY_UNIT)
    num_dens = C * AVOGADRO # mol/L * N/mol = N/L
    return uconvert(DEFAULT_NUMBERDENSITY_UNIT, num_dens)
end


# --- Molality Conversions ---
function cconvert(
    m_val, ::Type{Molality}, ::Type{Molarity};
    M_solute, rho_solution, kwargs...
)
    m = _ensure_unit(m_val, DEFAULT_MOLALITY_UNIT)
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)
    rho_s = _ensure_unit(rho_solution, DEFAULT_DENSITY_UNIT)

    m_solvent_basis = 1u"kg" # Basis mass
    n_solute = uconvert(u"mol", m * m_solvent_basis)
    m_solute = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solute * Ms)
    m_solution = m_solvent_basis + m_solute

    if m_solution <= 0 * unit(m_solution) throw(ArgumentError("Non-positive solution mass ($m_solution)")) end
    V_solution = uconvert(DEFAULT_VOLUME_UNIT, m_solution / rho_s)
    if V_solution <= 0 * unit(V_solution) throw(ArgumentError("Non-positive solution volume ($V_solution)")) end

    molarity = n_solute / V_solution
    return uconvert(DEFAULT_MOLARITY_UNIT, molarity)
end

function cconvert(
    m_val, ::Type{Molality}, ::Type{MoleFraction};
    M_solvent, kwargs...
)
    m = _ensure_unit(m_val, DEFAULT_MOLALITY_UNIT)
    Msv = _ensure_unit(M_solvent, DEFAULT_MASS_UNIT)

    m_solvent_basis = 1u"kg" # Basis mass
    n_solute = uconvert(u"mol", m * m_solvent_basis)
    n_solvent = uconvert(u"mol", m_solvent_basis / Msv)
    n_total = n_solute + n_solvent

    return n_total == 0 * unit(n_total) ? 0.0 : ustrip(NoUnits, n_solute / n_total)
end

function cconvert(
    m_val, ::Type{Molality}, ::Type{MassFraction};
    M_solute, kwargs...
)
    m = _ensure_unit(m_val, DEFAULT_MOLALITY_UNIT)
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)

    m_solvent_basis = 1u"kg" # Basis mass
    n_solute = uconvert(u"mol", m * m_solvent_basis)
    m_solute = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solute * Ms)
    m_solution = m_solvent_basis + m_solute

    if m_solution <= 0 * unit(m_solution) throw(ArgumentError("Non-positive solution mass ($m_solution)")) end

    mf = m_solute / m_solution
    return ustrip(NoUnits, mf)
end

# Molality -> VolumeFraction and Molality -> NumberDensity require rho_solution
function cconvert(
    m_val, ::Type{Molality}, ::Type{VolumeFraction};
    M_solute, rho_solute, rho_solution, kwargs...
)
    # Need to find V_solute and V_solution for the same amount of substance.
    # Go via Molarity first?
    C = cconvert(m_val, Molality, Molarity; M_solute, rho_solution)
    return cconvert(C, Molarity, VolumeFraction; M_solute, rho_solute)
end

function cconvert(
    m_val, ::Type{Molality}, ::Type{NumberDensity};
    M_solute, rho_solution, kwargs...
)
    # Go via Molarity first
    C = cconvert(m_val, Molality, Molarity; M_solute, rho_solution)
    return cconvert(C, Molarity, NumberDensity)
end

# --- Mole Fraction Conversions ---
function cconvert(
    chi_val::Real, ::Type{MoleFraction}, ::Type{Molarity};
    M_solute, M_solvent, rho_solution, kwargs...
)
    if !(0 <= chi_val <= 1) throw(ArgumentError("Mole fraction must be 0-1")) end
    chi = chi_val * NoUnits

    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)
    Msv = _ensure_unit(M_solvent, DEFAULT_MASS_UNIT)
    rho_s = _ensure_unit(rho_solution, DEFAULT_DENSITY_UNIT)

    n_total_basis = 1u"mol"
    n_solute = chi * n_total_basis
    n_solvent = (1*NoUnits - chi) * n_total_basis

    m_solute = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solute * Ms)
    m_solvent = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solvent * Msv)
    m_solution = m_solute + m_solvent

    if m_solution <= 0 * unit(m_solution) throw(ArgumentError("Non-positive solution mass ($m_solution)")) end
    V_solution = uconvert(DEFAULT_VOLUME_UNIT, m_solution / rho_s)
    if V_solution <= 0 * unit(V_solution) throw(ArgumentError("Non-positive solution volume ($V_solution)")) end

    molarity = n_solute / V_solution
    return uconvert(DEFAULT_MOLARITY_UNIT, molarity)
end

function cconvert(
    chi_val::Real, ::Type{MoleFraction}, ::Type{Molality};
    M_solvent, kwargs...
)
    if !(0 <= chi_val < 1) throw(ArgumentError("Mole fraction must be 0 <= χ < 1 for Molality")) end
    chi = chi_val * NoUnits
    Msv = _ensure_unit(M_solvent, DEFAULT_MASS_UNIT)

    n_total_basis = 1u"mol"
    n_solute = chi * n_total_basis
    n_solvent = (1*NoUnits - chi) * n_total_basis

    m_solvent = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solvent * Msv)
    if m_solvent <= 0 * unit(m_solvent) throw(ArgumentError("Non-positive solvent mass ($m_solvent), χ=1?")) end

    molality = n_solute / m_solvent
    return uconvert(DEFAULT_MOLALITY_UNIT, molality)
end

function cconvert(
    chi_val::Real, ::Type{MoleFraction}, ::Type{MassFraction};
    M_solute, M_solvent, kwargs...
)
    if !(0 <= chi_val <= 1) throw(ArgumentError("Mole fraction must be 0-1")) end
    chi = chi_val * NoUnits
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)
    Msv = _ensure_unit(M_solvent, DEFAULT_MASS_UNIT)

    n_total_basis = 1u"mol"
    n_solute = chi * n_total_basis
    n_solvent = (1*NoUnits - chi) * n_total_basis

    m_solute = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solute * Ms)
    m_solvent = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solvent * Msv)
    m_solution = m_solute + m_solvent

    if m_solution <= 0 * unit(m_solution) throw(ArgumentError("Non-positive solution mass ($m_solution)")) end

    mf = m_solute / m_solution
    return ustrip(NoUnits, mf)
end

# MoleFraction -> VolumeFraction/NumberDensity require rho_solution
function cconvert(
    chi_val::Real, ::Type{MoleFraction}, ::Type{VolumeFraction};
    M_solute, M_solvent, rho_solute, rho_solution, kwargs...
)
    # Go via Molarity
    C = cconvert(chi_val, MoleFraction, Molarity; M_solute, M_solvent, rho_solution)
    return cconvert(C, Molarity, VolumeFraction; M_solute, rho_solute)
end

function cconvert(
    chi_val::Real, ::Type{MoleFraction}, ::Type{NumberDensity};
    M_solute, M_solvent, rho_solution, kwargs...
)
    # Go via Molarity
    C = cconvert(chi_val, MoleFraction, Molarity; M_solute, M_solvent, rho_solution)
    return cconvert(C, Molarity, NumberDensity)
end

# --- Mass Fraction Conversions ---
function cconvert(
    mf_val::Real, ::Type{MassFraction}, ::Type{Molarity};
    M_solute, rho_solution, kwargs...
)
    if !(0 <= mf_val <= 1) throw(ArgumentError("Mass fraction must be 0-1")) end
    mf = mf_val * NoUnits
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)
    rho_s = _ensure_unit(rho_solution, DEFAULT_DENSITY_UNIT)

    m_solution_basis = 1u"kg" # Basis mass
    m_solute = mf * m_solution_basis
    n_solute = uconvert(u"mol", m_solute / Ms)

    V_solution = uconvert(DEFAULT_VOLUME_UNIT, m_solution_basis / rho_s)
    if V_solution <= 0 * unit(V_solution) throw(ArgumentError("Non-positive solution volume ($V_solution)")) end

    molarity = n_solute / V_solution
    return uconvert(DEFAULT_MOLARITY_UNIT, molarity)
end

function cconvert(
    mf_val::Real, ::Type{MassFraction}, ::Type{Molality};
    M_solute, kwargs...
)
    if !(0 <= mf_val < 1) throw(ArgumentError("Mass fraction must be 0 <= mf < 1 for Molality")) end
    mf = mf_val * NoUnits
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)

    m_solution_basis = 1u"kg" # Basis mass
    m_solute = mf * m_solution_basis
    m_solvent = (1*NoUnits - mf) * m_solution_basis

    if m_solvent <= 0 * unit(m_solvent) throw(ArgumentError("Non-positive solvent mass ($m_solvent), mf=1?")) end
    n_solute = uconvert(u"mol", m_solute / Ms)

    molality = n_solute / m_solvent
    return uconvert(DEFAULT_MOLALITY_UNIT, molality)
end

function cconvert(
    mf_val::Real, ::Type{MassFraction}, ::Type{MoleFraction};
    M_solute, M_solvent, kwargs...
)
    if !(0 <= mf_val <= 1) throw(ArgumentError("Mass fraction must be 0-1")) end
    mf = mf_val * NoUnits
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)
    Msv = _ensure_unit(M_solvent, DEFAULT_MASS_UNIT)

    m_solution_basis = 1u"kg" # Basis mass
    m_solute = mf * m_solution_basis
    m_solvent = (1*NoUnits - mf) * m_solution_basis

    n_solute = uconvert(u"mol", m_solute / Ms)
    # Handle mf=1 case (pure solute) -> n_solvent = 0
    n_solvent = if m_solvent > 0*unit(m_solvent) uconvert(u"mol", m_solvent / Msv) else 0u"mol" end
    n_total = n_solute + n_solvent

    return n_total == 0u"mol" ? 0.0 : ustrip(NoUnits, n_solute / n_total)
end

function cconvert(
    mf_val::Real, ::Type{MassFraction}, ::Type{VolumeFraction};
    rho_solute, rho_solution, kwargs...
)
    # Note: This uses rho_solution, relating V_solute to V_solution via masses
    if !(0 <= mf_val <= 1) throw(ArgumentError("Mass fraction must be 0-1")) end
    mf = mf_val * NoUnits
    rho_pure_s = _ensure_unit(rho_solute, DEFAULT_DENSITY_UNIT)
    rho_s = _ensure_unit(rho_solution, DEFAULT_DENSITY_UNIT)

    m_solution_basis = 1u"kg" # Basis mass
    m_solute = mf * m_solution_basis

    V_solute = uconvert(DEFAULT_VOLUME_UNIT, m_solute / rho_pure_s)
    V_solution = uconvert(DEFAULT_VOLUME_UNIT, m_solution_basis / rho_s)
    if V_solution <= 0 * unit(V_solution) throw(ArgumentError("Non-positive solution volume ($V_solution)")) end

    vf = V_solute / V_solution
    return ustrip(NoUnits, vf)
end

function cconvert(
    mf_val::Real, ::Type{MassFraction}, ::Type{NumberDensity};
    M_solute, rho_solution, kwargs...
)
    # Go via Molarity
    C = cconvert(mf_val, MassFraction, Molarity; M_solute, rho_solution)
    return cconvert(C, Molarity, NumberDensity)
end

# --- Volume Fraction Conversions ---
function cconvert(
    vf_val::Real, ::Type{VolumeFraction}, ::Type{Molarity};
    M_solute, rho_solute, kwargs...
)
    # Inverse of Molarity -> VolumeFraction
    # Assumes vf = V_solute / V_solution
    if !(0 <= vf_val <= 1) throw(ArgumentError("Volume fraction must be 0-1")) end
    vf = vf_val * NoUnits
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)
    rho_pure_s = _ensure_unit(rho_solute, DEFAULT_DENSITY_UNIT)

    V0 = 1u"L" # Basis volume
    V_solute = vf * V0
    m_solute = uconvert(DEFAULT_MASS_QUANTITY_UNIT, V_solute * rho_pure_s)
    n_solute = uconvert(u"mol", m_solute / Ms)

    molarity = n_solute / V0
    return uconvert(DEFAULT_MOLARITY_UNIT, molarity)
end

function cconvert(
    vf_val::Real, ::Type{VolumeFraction}, ::Type{MassFraction};
    rho_solute, rho_solution, kwargs...
)
    # Inverse of MassFraction -> VolumeFraction
     if !(0 <= vf_val <= 1) throw(ArgumentError("Volume fraction must be 0-1")) end
     vf = vf_val * NoUnits
     rho_pure_s = _ensure_unit(rho_solute, DEFAULT_DENSITY_UNIT)
     rho_s = _ensure_unit(rho_solution, DEFAULT_DENSITY_UNIT)

     V_solution_basis = 1u"L" # Basis volume
     V_solute = vf * V_solution_basis

     m_solute = uconvert(DEFAULT_MASS_QUANTITY_UNIT, V_solute * rho_pure_s)
     m_solution = uconvert(DEFAULT_MASS_QUANTITY_UNIT, V_solution_basis * rho_s)
     if m_solution <= 0 * unit(m_solution) throw(ArgumentError("Non-positive solution mass ($m_solution)")) end

     mf = m_solute / m_solution
     return ustrip(NoUnits, mf)
end

# VolumeFraction -> Molality/MoleFraction/NumberDensity require more info or intermediate steps
function cconvert(
    vf_val::Real, ::Type{VolumeFraction}, ::Type{Molality};
    M_solute, rho_solute, rho_solution, kwargs...
)
    # Go via Molarity
    C = cconvert(vf_val, VolumeFraction, Molarity; M_solute, rho_solute)
    # Need M_solute and rho_solution for Molarity -> Molality step (already required)
    return cconvert(C, Molarity, Molality; M_solute, rho_solution)
end

function cconvert(
    vf_val::Real, ::Type{VolumeFraction}, ::Type{MoleFraction};
    M_solute, M_solvent, rho_solute, rho_solution, kwargs...
)
    # Go via Molarity
    C = cconvert(vf_val, VolumeFraction, Molarity; M_solute, rho_solute)
    # Need M_solute, M_solvent, rho_solution for Molarity -> MoleFraction step (already required)
    return cconvert(C, Molarity, MoleFraction; M_solute, M_solvent, rho_solution)
end

function cconvert(
    vf_val::Real, ::Type{VolumeFraction}, ::Type{NumberDensity};
    M_solute, rho_solute, kwargs...
)
     # Go via Molarity
    C = cconvert(vf_val, VolumeFraction, Molarity; M_solute, rho_solute)
    return cconvert(C, Molarity, NumberDensity)
end

# --- Number Density Conversions ---
function cconvert(N_val, ::Type{NumberDensity}, ::Type{Molarity}; kwargs...)
    N = _ensure_unit(N_val, DEFAULT_NUMBERDENSITY_UNIT)
    molarity = N / AVOGADRO
    return uconvert(DEFAULT_MOLARITY_UNIT, molarity)
end

# Other NumberDensity conversions go via Molarity
function cconvert(
    N_val, ::Type{NumberDensity}, ::Type{Molality};
    M_solute, rho_solution, kwargs...
)
    C = cconvert(N_val, NumberDensity, Molarity)
    return cconvert(C, Molarity, Molality; M_solute, rho_solution)
end

function cconvert(
    N_val, ::Type{NumberDensity}, ::Type{MoleFraction};
    M_solute, M_solvent, rho_solution, kwargs...
)
    C = cconvert(N_val, NumberDensity, Molarity)
    return cconvert(C, Molarity, MoleFraction; M_solute, M_solvent, rho_solution)
end

function cconvert(
    N_val, ::Type{NumberDensity}, ::Type{MassFraction};
    M_solute, rho_solution, kwargs...
)
    C = cconvert(N_val, NumberDensity, Molarity)
    return cconvert(C, Molarity, MassFraction; M_solute, rho_solution)
end

function cconvert(
    N_val, ::Type{NumberDensity}, ::Type{VolumeFraction};
    M_solute, rho_solute, kwargs...
)
    C = cconvert(N_val, NumberDensity, Molarity)
    return cconvert(C, Molarity, VolumeFraction; M_solute, rho_solute)
end

const UNIT_TYPE_MAP = Dict(
    u"mol/L" => "mol/L",
    u"mol/kg" => "mol/kg",
    u"1/(Å^3)" => "molecule/Å^3",
)

# ==============================================================================
# User-Friendly Wrapper (String-based, handles percentages)
# ==============================================================================
public UNIT_STRINGS

"""
    UNIT_STRINGS

A mapping of string identifiers to concentration unit types.

```julia
julia> Packmol.UNIT_STRINGS
OrderedCollections.OrderedDict{String, DataType} with 26 entries:
  "Molarity"       => Molarity
  "molarity"       => Molarity
  "mol/L"          => Molarity
  "molar"          => Molarity
  "M"              => Molarity
  "Molality"       => Molality
  "molality"       => Molality
  "mol/kg"         => Molality
  "molal"          => Molality
  "molefraction"   => MoleFraction
  "chi"            => MoleFraction
  "x"              => MoleFraction
  "MoleFraction"   => MoleFraction
  "mole fraction"  => MoleFraction
  "MassFraction"   => MassFraction
  "massfraction"   => MassFraction
  "%m/m"           => MassFraction
  "w/w"            => MassFraction
  "mass %"         => MassFraction
  "%w/w"           => MassFraction
  "VolumeFraction" => VolumeFraction
  "volumefraction" => VolumeFraction
  "%v/v"           => VolumeFraction
  "v/v"            => VolumeFraction
  "volume %"       => VolumeFraction
  "%vol/vol"       => VolumeFraction
  "NumberDensity"  => NumberDensity
  "numberdensity"  => NumberDensity
  "Å^-3"           => NumberDensity,
  "molecule/angs3" => NumberDensity,
  "molecules/angs3"=> NumberDensity,
  "molecule/Å^-3"  => NumberDensity,
  "molecules/Å^-3" => NumberDensity,
  "molecule/Å³"    => NumberDensity,
  "molecules/Å³"   => NumberDensity,
```

"""
UNIT_STRINGS

const UNIT_STRINGS = OrderedDict(
  "Molarity"       => Molarity,
  "molarity"       => Molarity,
  "mol/L"          => Molarity,
  "molar"          => Molarity,
  "M"              => Molarity,
  "Molality"       => Molality,
  "molality"       => Molality,
  "mol/kg"         => Molality,
  "molal"          => Molality,
  "molefraction"   => MoleFraction,
  "chi"            => MoleFraction,
  "x"              => MoleFraction,
  "MoleFraction"   => MoleFraction,
  "mole fraction"  => MoleFraction,
  "MassFraction"   => MassFraction,
  "massfraction"   => MassFraction,
  "%m/m"           => MassFraction,
  "w/w"            => MassFraction,
  "mass %"         => MassFraction,
  "%w/w"           => MassFraction,
  "VolumeFraction" => VolumeFraction,
  "volumefraction" => VolumeFraction,
  "%v/v"           => VolumeFraction,
  "v/v"            => VolumeFraction,
  "volume %"       => VolumeFraction,
  "%vol/vol"       => VolumeFraction,
  "NumberDensity"  => NumberDensity,
  "numberdensity"  => NumberDensity,
  "Å^-3"           => NumberDensity,
  "molecule/angs3" => NumberDensity,
  "molecules/angs3"=> NumberDensity,
  "molecule/Å^-3"  => NumberDensity,
  "molecules/Å^-3" => NumberDensity,
  "molecule/Å³"    => NumberDensity,
  "molecules/Å³"   => NumberDensity,
)

"""
    cconvert(value, units::Pair{String, String}; kwargs...)

Convert a concentration `value` from one unit to another using string identifiers.
Handles conversion between percentage and fraction representations for dimensionless units.

# Arguments
- `value`: The numerical value (or `Unitful.Quantity`) of the input concentration.
- `units`: A `Pair` of strings specifying the conversion, e.g., `"mol/L" => "mol/kg"` or `"10 %m/m" => "MoleFraction"`.
           Supported unit strings include common abbreviations and names (see `Packmol.UNIT_STRINGS`).
           For dimensionless units (`%m/m`, `%v/v`), providing a "%" sign implies the input `value` is a percentage (e.g., 10.0 for 10%).
           Otherwise, it's treated as a fraction (0-1). Output format also depends on "%" in the target string.
- `kwargs`: Keyword arguments providing necessary auxiliary data (e.g., `M_solute`, `rho_solution`).
            These can be numbers (assuming default units like g/mol, kg/L) or `Unitful.Quantity` objects.

# Returns
- The converted concentration. This will be a `Unitful.Quantity` for Molarity, Molality, NumberDensity,
  or a `Real` number for dimensionless fractions (or percentages if requested).

# Examples

```julia-repl
julia> using Packmol

# Example: Ethanol (EtOH) in Water (H2O) mixture
M_EtOH = 46.068u"g/mol"
M_H2O = 18.015u"g/mol"
rho_EtOH_pure = 0.789u"kg/L" # Density of pure Ethanol
# Some realistic solution densities for examples:
rho_sol_10M = 0.90u"kg/L"   # Approx. density for ~10 mol/L EtOH (~58% w/w)
rho_sol_50ww = 0.914u"kg/L" # Approx. density for 50% w/w EtOH
rho_sol_chi02 = 0.94u"kg/L"  # Approx. density for mole fraction χ_EtOH = 0.2

println("--- Ethanol/Water Examples ---")

# Molarity to Molality (~10 M solution)
C_in1 = 10.0u"mol/L"
m1 = cconvert(C_in1, "mol/L" => "mol/kg"; M_solute=M_EtOH, rho_solution=rho_sol_10M) # Expect ~22.76 mol/kg

# Mass Percent (input %) to Molarity (output Quantity)
mass_perc_in = 50.0
C1 = cconvert(mass_perc_in, "%m/m" => "Molarity"; M_solute=M_EtOH, rho_solution=rho_sol_50ww) # Expect ~9.92 mol/L

# Mole fraction (input fraction) to Molality (output Quantity)
chi_in = 0.2
m2 = cconvert(chi_in, "mole fraction" => "molality"; M_solvent=M_H2O) # Expect ~13.88 mol/kg

# Molarity (input Quantity) to Volume Percent (output %)
C_in2 = 10.0u"M" # Unitful recognizes M as mol/L
vp1 = cconvert(C_in2, "M" => "%v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) # Expect ~58.39 %v/v

# Number Density (input Quantity) to Molarity
N_in = 0.005u"Å^-3"
C2 = cconvert(N_in, "Å^-3" => "Molarity") # Expect ~8.30 mol/L

# Example using plain numbers for kwargs (assuming default units)
m3 = cconvert(0.2, "mole fraction" => "molality"; M_solvent=18.015) # M_H2O in g/mol

```

"""
function cconvert(value, units::Pair{String, String}; kwargs...)
    unit_in_str = strip(units.first)
    unit_out_str = strip(units.second)

    T_in = get(UNIT_STRINGS, unit_in_str, nothing)
    T_out = get(UNIT_STRINGS, unit_out_str, nothing)
    
    if value isa Quantity
        oneunit(value) == _oneunit(T_in) || throw(ArgumentError("Unit of input value ($value) does not match input unit type $T_in ($unit_in_str)"))
    end
    isnothing(T_in) && throw(ArgumentError("Unknown input concentration unit: $(units.first)"))
    isnothing(T_out) && throw(ArgumentError("Unknown output concentration unit: $(units.second)"))

    # Handle % vs fraction for dimensionless units based on input string
    input_val = value
    is_percent_in = occursin('%', unit_in_str) || occursin("percent", unit_in_str) # Check original case for %? No, lower is fine.
    is_fraction_type_in = T_in in (MassFraction, VolumeFraction, MoleFraction)
    
    if is_percent_in && is_fraction_type_in
        if input_val isa Quantity # If user passed 10.0u"percent" or similar
            input_val = ustrip(input_val) / 100.0
        elseif input_val isa Real
            input_val = input_val / 100.0 # Assume input 10.0 means 10% -> 0.1 fraction
        end
    end
    
    # Call the type-dispatch version - ensure kwargs are passed
    # Need to handle potential UndefKeywordErrors gracefully if required kwargs are missing.
    result = cconvert(input_val, T_in, T_out; kwargs...)
    
    # Handle % vs fraction for output based on output string
    is_percent_out = occursin('%', unit_out_str) || occursin("percent", unit_out_str)
    is_fraction_type_out = T_out in (MassFraction, VolumeFraction, MoleFraction)
    
    if is_percent_out && is_fraction_type_out
        # Result should be a Real number (dimensionless fraction) from internal cconvert
        if result isa Quantity && dimension(result) == NoDims
            return ustrip(result) * 100.0
        elseif result isa Real
            return result * 100.0
        else
            @error "Internal conversion to fraction type $T_out returned unexpected type $(typeof(result)). Cannot convert to percentage."
            return result # Return the unexpected result
        end
    else
        # Return raw fraction (Real) or Quantity with units
        return result
    end

end

@testsnippet CConvert begin
    using Packmol
    using Unitful: ustrip, Na, dimensiohn, Quantity, uconvert

    # --- Test Data: Ethanol (EtOH) in Water (H2O) mixture at 25oC
    M_EtOH = 46.068u"g/mol"
    M_H2O = 18.015u"g/mol"
    rho_EtOH_pure = 0.785u"kg/L"
    rho_H2O_pure = 0.997u"kg/L"
    
    # Realistic solution densities for specific concentrations at 25oC
    rho_sol_10M = 0.845u"kg/L"    # Approx. density for ~10 mol/L EtOH
    rho_sol_50ww = 0.914u"kg/L"  # Approx. density for 50% w/w EtOH
    rho_sol_50vv = 0.93u"kg/L"   # Approx. density for 50% v/v EtOH
    rho_sol_chi02 = 0.932u"kg/L"   # Approx. density for mole fraction χ_EtOH = 0.2
    rho_sol_1M = 0.982u"kg/L"   # Approx. density for 1 mol/L EtOH (~4.5% w/w)
    rho_sol_1m = 0.992u"kg/L"   # Approx. density for 1 mol/kg EtOH (~4.4% w/w)

    #voltar
end

@testitem "Molarity Conversions" setup=[CConvert] begin

    # --- Identity ---
    @test cconvert(1.0u"mol/L", "Molarity" => "Molarity") == 1.0u"mol/L"
    @test cconvert(1.0, "M" => "mol/L") == 1.0u"mol/L" # Number input

    # --- Molarity to Molality ---
    # Simple case: C=1, Ms=100, rho=1 -> m = 1/(1-0.1) = 1/0.9 = 1.111...
    @test cconvert(1.0u"mol/L", "M" => "mol/kg"; M_solute=100u"g/mol", rho_solution=1.0u"kg/L") ≈ 1.11u"mol/kg" rtol=1e-3
    # Ethanol case (~10 M -> ~22.76 m)
    @test cconvert(10.0u"M", "M" => "mol/kg"; M_solute=M_EtOH, rho_solution=rho_sol_10M) ≈ 26.019u"mol/kg" rtol=1e-3
    # Ethanol case (1 M -> ~1.036 m)
    @test cconvert(1.0u"M", "M" => "mol/kg"; M_solute=M_EtOH, rho_solution=rho_sol_1M) ≈ 1.0684u"mol/kg" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "mol/kg"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "mol/kg"; rho_solution=rho_sol_1M)
    # Zero concentration
    @test cconvert(0.0u"M", "M" => "mol/kg"; M_solute=M_EtOH, rho_solution=rho_sol_1M) == 0.0u"mol/kg"

    # --- Molarity to MoleFraction ---
    # Ethanol case (1 M -> ~0.0183)
    @test cconvert(1.0u"M", "M" => "MoleFraction"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_1M) ≈ 0.01888 rtol=1e-3
    # Ethanol case (~10 M -> ~0.359)
    @test cconvert(10.0u"M", "M" => "MoleFraction"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_10M) ≈ 0.31914 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O)
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "chi"; M_solute=M_EtOH, rho_solution=rho_sol_1M)
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "chi"; M_solvent=M_H2O, rho_solution=rho_sol_1M)
    # Zero concentration
    @test cconvert(0.0u"M", "M" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_1M) == 0.0

    # --- Molarity to MassFraction ---
    # Ethanol case (1 M -> ~0.04686)
    @test cconvert(1.0u"M", "M" => "MassFraction"; M_solute=M_EtOH, rho_solution=rho_sol_1M) ≈ (1.0*ustrip(M_EtOH)/1000)/ustrip(rho_sol_1M) # 0.046068 / 0.983
    @test cconvert(1.0u"M", "M" => "MassFraction"; M_solute=M_EtOH, rho_solution=rho_sol_1M) ≈ 0.04691 rtol=1e-2
    # Ethanol case (~10 M -> ~0.5118)
    @test cconvert(10.0u"M", "M" => "w/w"; M_solute=M_EtOH, rho_solution=rho_sol_10M) ≈ (10.0*ustrip(M_EtOH)/1000)/ustrip(rho_sol_10M) # 0.46068 / 0.90
    @test cconvert(10.0u"M", "M" => "w/w"; M_solute=M_EtOH, rho_solution=rho_sol_10M) ≈ 0.54518 rtol=1e-2
    # Output as Percentage
    @test cconvert(1.0u"M", "M" => "%m/m"; M_solute=M_EtOH, rho_solution=rho_sol_1M) ≈ 4.69124 rtol=1e-2
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "w/w"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "w/w"; rho_solution=rho_sol_1M)
    # Zero concentration
    @test cconvert(0.0u"M", "M" => "w/w"; M_solute=M_EtOH, rho_solution=rho_sol_1M) == 0.0
    @test cconvert(0.0u"M", "M" => "%m/m"; M_solute=M_EtOH, rho_solution=rho_sol_1M) == 0.0

    # --- Molarity to VolumeFraction ---
    # Ethanol case (1 M -> ~0.05839)
    @test cconvert(1.0u"M", "M" => "VolumeFraction"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ (1.0*ustrip(M_EtOH)/1000)/ustrip(rho_EtOH_pure) # 0.046068/0.789
    @test cconvert(1.0u"M", "M" => "VolumeFraction"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 0.05868 rtol=1e-3
    # Ethanol case (~10 M -> ~0.5839)
    @test cconvert(10.0u"M", "M" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 0.5868 rtol=1e-3
    # Output as Percentage
    @test cconvert(10.0u"M", "M" => "%v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 58.68 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "v/v"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(1.0u"M", "M" => "v/v"; rho_solute=rho_EtOH_pure)
    # Zero concentration
    @test cconvert(0.0u"M", "M" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) == 0.0
    @test cconvert(0.0u"M", "M" => "%v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) == 0.0

    # --- Molarity to NumberDensity ---
    @test cconvert(1.0u"mol/L", "M" => "Å^-3") ≈ 1.0u"mol/L" * Na rtol=1e-6
    @test cconvert(1.0u"mol/L", "M" => "Å^-3") ≈ 0.000602214u"Å^-3"  rtol=1e-6 
    # Missing Kwargs (none needed)
    @test cconvert(1.0u"M", "M" => "Å^-3") isa Quantity
    # Zero concentration
    @test cconvert(0.0u"M", "M" => "Å^-3") == 0.0u"Å^-3"

end

@testitem "Molality Conversions" setup=[CConvert] begin

    # --- Identity ---
    @test cconvert(1.0u"mol/kg", "Molality" => "Molality") == 1.0u"mol/kg"
    @test cconvert(1.0, "mol/kg" => "molal") == 1.0u"mol/kg" # Number input

    # --- Molality to Molarity ---
    # Ethanol case (1 m -> ~0.97 M) - use density for 1m solution
    @test cconvert(1.0u"mol/kg", "mol/kg" => "M"; M_solute=M_EtOH, rho_solution=rho_sol_1m) ≈ 0.9483u"mol/L" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "M"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "M"; rho_solution=rho_sol_1m)
    # Zero concentration
    @test cconvert(0.0u"mol/kg", "mol/kg" => "M"; M_solute=M_EtOH, rho_solution=rho_H2O_pure) == 0.0u"mol/L" # Use water density for zero conc

    # --- Molality to MoleFraction ---
    # Ethanol case (1 m -> ~0.0177)
    @test cconvert(1.0u"mol/kg", "mol/kg" => "chi"; M_solvent=M_H2O) ≈ 1.0 / (1.0 + 1000/ustrip(M_H2O))
    @test cconvert(1.0u"mol/kg", "mol/kg" => "chi"; M_solvent=M_H2O) ≈ 0.017696 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "chi")
    # Zero concentration
    @test cconvert(0.0u"mol/kg", "mol/kg" => "chi"; M_solvent=M_H2O) == 0.0

    # --- Molality to MassFraction ---
    # Ethanol case (1 m -> ~0.0440)
    m_s = 1.0 * ustrip(M_EtOH)
    mf = m_s / (1000 + m_s)
    @test cconvert(1.0u"mol/kg", "mol/kg" => "w/w"; M_solute=M_EtOH) ≈ mf
    @test cconvert(1.0u"mol/kg", "mol/kg" => "w/w"; M_solute=M_EtOH) ≈ 0.044047 rtol=1e-3
    # Output Percentage
    @test cconvert(1.0u"mol/kg", "mol/kg" => "%m/m"; M_solute=M_EtOH) ≈ 4.4047 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "w/w")
    # Zero concentration
    @test cconvert(0.0u"mol/kg", "mol/kg" => "w/w"; M_solute=M_EtOH) == 0.0
    @test cconvert(0.0u"mol/kg", "mol/kg" => "%m/m"; M_solute=M_EtOH) == 0.0

    # --- Molality to VolumeFraction --- (Requires rho_solution, rho_solute)
    @test cconvert(1.0u"mol/kg", "mol/kg" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_1m) ≈ 0.0557 rtol=1e-3
    # Output Percentage
    @test cconvert(1.0u"mol/kg", "mol/kg" => "%v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_1m) ≈ 5.5652 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) # Missing rho_solution
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "v/v"; M_solute=M_EtOH, rho_solution=rho_sol_1m) # Missing rho_solute
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_1m) # Missing M_solute
    # Zero concentration
    @test cconvert(0.0u"mol/kg", "mol/kg" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure, rho_solution=rho_H2O_pure) == 0.0

    # --- Molality to NumberDensity --- (Requires rho_solution)
    # Ethanol case (1m -> intermediate molarity -> number density)
    nd = cconvert(1.0u"mol/kg", "mol/kg" => "Å^-3"; M_solute=M_EtOH, rho_solution=rho_sol_1m)
    C_interim = cconvert(1.0u"mol/kg", "mol/kg" => "M"; M_solute=M_EtOH, rho_solution=rho_sol_1m)
    @test nd ≈ cconvert(C_interim, "mol/L" => "Å^-3")
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "Å^-3"; M_solute=M_EtOH) # Missing rho_solution
    @test_throws UndefKeywordError cconvert(1.0u"mol/kg", "mol/kg" => "Å^-3"; rho_solution=rho_sol_1m) # Missing M_solute
    # Zero concentration
    @test cconvert(0.0u"mol/kg", "mol/kg" => "Å^-3"; M_solute=M_EtOH, rho_solution=rho_H2O_pure) == 0.0u"Å^-3"

end

@testitem "MoleFraction Conversions" setup=[CConvert] begin

    # --- Identity ---
    @test cconvert(0.5, "MoleFraction" => "MoleFraction") == 0.5
    @test cconvert(0.5, "chi" => "x") == 0.5 # String variations

    # --- MoleFraction to Molarity ---
    # Ethanol case (χ=0.2 -> ~5.71 M) - use density for χ=0.2 solution
    @test cconvert(0.2, "chi" => "M"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_chi02) ≈ 7.88974u"mol/L" rtol=1e-3
    @test cconvert(1.0, "chi" => "M"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_EtOH_pure) ≈ 17.040u"mol/L" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "M"; M_solute=M_EtOH, M_solvent=M_H2O)
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "M"; M_solute=M_EtOH, rho_solution=rho_sol_chi02)
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "M"; M_solvent=M_H2O, rho_solution=rho_sol_chi02)
    # Edge Cases
    @test cconvert(0.0, "chi" => "M"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_H2O_pure) == 0.0u"mol/L"

    # --- MoleFraction to Molality ---
    # Ethanol case (χ=0.2 -> ~13.88 m)
    @test cconvert(0.2, "chi" => "mol/kg"; M_solvent=M_H2O) ≈ (0.2 / (0.8 * ustrip(M_H2O)/1000.0)) * u"mol/kg"
    @test cconvert(0.2, "chi" => "mol/kg"; M_solvent=M_H2O) ≈ 13.877u"mol/kg" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "mol/kg")
    # Edge Cases
    @test cconvert(0.0, "chi" => "mol/kg"; M_solvent=M_H2O) == 0.0u"mol/kg"
    @test_throws ArgumentError cconvert(1.0, "chi" => "mol/kg"; M_solvent=M_H2O) # Requires chi < 1

    # --- MoleFraction to MassFraction ---
    # Ethanol case (χ=0.2 -> ~0.390)
    m_etoh = 0.2 * ustrip(M_EtOH)
    m_h2o = 0.8 * ustrip(M_H2O)
    mf = m_etoh / (m_etoh + m_h2o)
    @test cconvert(0.2, "chi" => "w/w"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ mf
    @test cconvert(0.2, "chi" => "w/w"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ 0.38998 rtol=1e-3
    # Output Percentage
    @test cconvert(0.2, "chi" => "%m/m"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ 38.998 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "w/w"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "w/w"; M_solvent=M_H2O)
    # Edge Cases
    @test cconvert(0.0, "chi" => "w/w"; M_solute=M_EtOH, M_solvent=M_H2O) == 0.0
    @test cconvert(1.0, "chi" => "w/w"; M_solute=M_EtOH, M_solvent=M_H2O) == 1.0
    @test cconvert(1.0, "chi" => "%m/m"; M_solute=M_EtOH, M_solvent=M_H2O) == 100.0

    # --- MoleFraction to VolumeFraction --- (Requires densities)
    # Ethanol case (χ=0.2 -> ~0.265) - Requires intermediate Molarity calc
    @test cconvert(0.2, "chi" => "v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_chi02) ≈ 0.4630 rtol=1e-3
    # Output Percentage
    @test cconvert(0.2, "chi" => "%v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_chi02) ≈ 46.30 rtol=1e-3 
    # Missing Kwargs (will cascade from intermediate call)
    @test_throws UndefKeywordError cconvert(0.2, "chi" => "v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure) # Missing rho_solution
    # Edge cases
    @test cconvert(0.0, "chi" => "v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_H2O_pure) == 0.0
    # Pure solute vf=1 needs density of pure solute as rho_solution
    @test cconvert(1.0, "chi" => "v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_EtOH_pure) ≈ 1.0

     # --- MoleFraction to NumberDensity --- (Requires densities)
     # Ethanol case (χ=0.2 -> ~3.44e-3 Å^-3) - Requires intermediate Molarity calc
     @test cconvert(0.2, "chi" => "Å^-3"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_chi02) ≈ 0.0047513u"Å^-3" rtol=1e-3
     # Missing Kwargs
     @test_throws UndefKeywordError cconvert(0.2, "chi" => "Å^-3"; M_solute=M_EtOH, M_solvent=M_H2O) # Missing rho_solution
     # Edge Cases
     @test cconvert(0.0, "chi" => "Å^-3"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_H2O_pure) == 0.0u"Å^-3"

end

@testitem "MassFraction Conversions" setup=[CConvert] begin

    # --- Identity ---
    @test cconvert(0.5, "MassFraction" => "MassFraction") == 0.5
    @test cconvert(50.0, "%m/m" => "w/w") == 0.5 # Percent input, fraction output
    @test cconvert(0.5, "w/w" => "%m/m") == 50.0 # Fraction input, percent output

    # --- MassFraction to Molarity ---
    # Ethanol case (50% w/w -> ~9.92 M)
    @test cconvert(50.0, "%m/m" => "M"; M_solute=M_EtOH, rho_solution=rho_sol_50ww) ≈ 9.920u"mol/L" rtol=1e-3
    @test cconvert(0.5, "w/w" => "M"; M_solute=M_EtOH, rho_solution=rho_sol_50ww) ≈ 9.920u"mol/L" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "M"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "M"; rho_solution=rho_sol_50ww)
    # Edge Cases
    @test cconvert(0.0, "w/w" => "M"; M_solute=M_EtOH, rho_solution=rho_H2O_pure) == 0.0u"mol/L"
    @test cconvert(100.0, "%m/m" => "M"; M_solute=M_EtOH, rho_solution=rho_EtOH_pure) ≈ uconvert(u"mol/L", rho_EtOH_pure / M_EtOH) # Molarity of pure ethanol
    @test cconvert(100.0, "%m/m" => "M"; M_solute=M_EtOH, rho_solution=rho_EtOH_pure) ≈ 17.040u"mol/L" rtol=1e-3

    # --- MassFraction to Molality ---
    # Ethanol case (50% w/w -> 21.71 m)
    @test cconvert(50.0, "%m/m" => "mol/kg"; M_solute=M_EtOH) ≈ (0.5 / ustrip(M_EtOH)) / 0.5 * 1000 * u"mol/kg" # (m_s/M_s)/m_sv
    @test cconvert(50.0, "%m/m" => "mol/kg"; M_solute=M_EtOH) ≈ 21.707u"mol/kg" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "mol/kg")
    # Edge Cases
    @test cconvert(0.0, "w/w" => "mol/kg"; M_solute=M_EtOH) == 0.0u"mol/kg"
    @test_throws ArgumentError cconvert(1.0, "w/w" => "mol/kg"; M_solute=M_EtOH) # mf must be < 1
    @test_throws ArgumentError cconvert(100.0, "%m/m" => "mol/kg"; M_solute=M_EtOH)

    # --- MassFraction to MoleFraction ---
    # Ethanol case (50% w/w -> χ ~ 0.259)
    n_etoh = 50 / ustrip(M_EtOH)
    n_h2o = 50 / ustrip(M_H2O)
    chi = n_etoh / (n_etoh + n_h2o)
    @test cconvert(50.0, "%m/m" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ chi
    @test cconvert(50.0, "%m/m" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ 0.28111 rtol=1e-3
    @test cconvert(50.0, "%m/m" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ 0.28111 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "chi"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "chi"; M_solvent=M_H2O)
    # Edge Cases
    @test cconvert(0.0, "w/w" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) == 0.0
    @test cconvert(1.0, "w/w" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) == 1.0
    @test cconvert(100.0, "%m/m" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) == 1.0

    # --- MassFraction to VolumeFraction ---
    # Ethanol case (50% w/w -> vf ~ 0.584)
    vf = (0.5 / ustrip(rho_EtOH_pure)) / (1.0 / ustrip(rho_sol_50ww))
    @test cconvert(50.0, "%m/m" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50ww) ≈ vf
    @test cconvert(50.0, "%m/m" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50ww) ≈ 0.5822 rtol=1e-3 
    @test cconvert(0.5, "w/w" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50ww) ≈ (0.5/ustrip(rho_EtOH_pure))/(1.0/ustrip(rho_sol_50ww))
    @test cconvert(0.5, "w/w" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50ww) ≈ 0.5822 rtol=1e-3
    # Output Percentage
    @test cconvert(50.0, "%m/m" => "%v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50ww) ≈ 58.216 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "v/v"; rho_solute=rho_EtOH_pure)
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "v/v"; rho_solution=rho_sol_50ww)
    # Edge Cases
    @test cconvert(0.0, "w/w" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_H2O_pure) == 0.0
    @test cconvert(1.0, "w/w" => "v/v"; rho_solute=rho_EtOH_pure, rho_solution=rho_EtOH_pure) == 1.0

    # --- MassFraction to NumberDensity ---
    # Ethanol case (50% w/w -> ~5.97e-3 Å^-3) - intermediate molarity
    @test cconvert(50.0, "%m/m" => "Å^-3"; M_solute=M_EtOH, rho_solution=rho_sol_50ww) ≈ 0.00597u"Å^-3"  rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "Å^-3"; M_solute=M_EtOH) # Missing rho_solution
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "Å^-3"; rho_solution=rho_sol_50ww) # Missing M_solute
    # Edge Cases
    @test cconvert(0.0, "w/w" => "Å^-3"; M_solute=M_EtOH, rho_solution=rho_H2O_pure) == 0.0u"Å^-3"

end

@testitem "VolumeFraction Conversions" setup=[CConvert] begin

    # --- Identity ---
    @test cconvert(0.5, "VolumeFraction" => "VolumeFraction") == 0.5
    @test cconvert(50.0, "%v/v" => "v/v") == 0.5 # Percent input, fraction output
    @test cconvert(0.5, "v/v" => "%v/v") == 50.0 # Fraction input, percent output

    # --- VolumeFraction to Molarity ---
    # Ethanol case (vf=0.5 -> ~8.56 M)
    @test cconvert(0.5, "v/v" => "M"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ (0.5 * ustrip(rho_EtOH_pure)) / (ustrip(M_EtOH)/1000.0) * u"mol/L"
    @test cconvert(0.5, "v/v" => "M"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 8.520u"mol/L" rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "M"; M_solute=M_EtOH)
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "M"; rho_solute=rho_EtOH_pure)
    # Edge Cases
    @test cconvert(0.0, "v/v" => "M"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) == 0.0u"mol/L"
    @test cconvert(1.0, "v/v" => "M"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ uconvert(u"mol/L", rho_EtOH_pure / M_EtOH) # Molarity of pure

    # --- VolumeFraction to MassFraction ---
    # Ethanol case (vf=0.5 -> mf ~ 0.424) - use rho for vf=0.5
    mf = (0.5 * ustrip(rho_EtOH_pure)) / ustrip(rho_sol_50vv)
    @test cconvert(0.5, "v/v" => "w/w"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50vv) ≈ mf
    @test cconvert(0.5, "v/v" => "w/w"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50vv) ≈ 0.4220 rtol=1e-3
    # Output Percentage
    @test cconvert(50.0, "%v/v" => "%m/m"; rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50vv) ≈ 42.20 rtol=1e-3
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "w/w"; rho_solute=rho_EtOH_pure)
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "w/w"; rho_solution=rho_sol_50vv)
    # Edge Cases
    @test cconvert(0.0, "v/v" => "w/w"; rho_solute=rho_EtOH_pure, rho_solution=rho_H2O_pure) == 0.0
    @test cconvert(1.0, "v/v" => "w/w"; rho_solute=rho_EtOH_pure, rho_solution=rho_EtOH_pure) == 1.0

    # --- VolumeFraction to Molality --- (Requires intermediate Molarity/MassFraction)
    # Ethanol case (vf=0.5 -> ~15.6 m)
    @test cconvert(0.5, "v/v" => "mol/kg"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50vv) ≈ 15.85u"mol/kg" rtol=1e-3 # Chain M->m or MF->m
    # Missing Kwargs (will cascade)
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "mol/kg"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) # Missing rho_solution
    # Edge Cases
    @test cconvert(0.0, "v/v" => "mol/kg"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure, rho_solution=rho_H2O_pure) == 0.0u"mol/kg"
    @test_throws ArgumentError cconvert(1.0, "v/v" => "mol/kg"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure, rho_solution=rho_EtOH_pure) # Fails at MF->m step (mf=1)

    # --- VolumeFraction to MoleFraction --- (Requires intermediate Molarity/MassFraction)
    # Ethanol case (vf=0.5 -> χ ~ 0.218)
    @test cconvert(0.5, "v/v" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_50vv) ≈ 0.2221 rtol=1e-3
    # Missing Kwargs (will cascade)
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure) # Missing rho_solution
    # Edge Cases
    @test cconvert(0.0, "v/v" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_H2O_pure) == 0.0
    @test cconvert(1.0 - eps(1.0), "v/v" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_EtOH_pure) ≈ 1.0

    # --- VolumeFraction to NumberDensity ---
    # Ethanol case (vf=0.5 -> ~5.16e-3 Å^-3) - Intermediate Molarity
    @test cconvert(0.5, "v/v" => "Å^-3"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 0.0051308u"Å^-3"  rtol=1e-4
    # Missing Kwargs
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "Å^-3"; M_solute=M_EtOH) # rho_solute missing
    @test_throws UndefKeywordError cconvert(0.5, "v/v" => "Å^-3"; rho_solute=rho_EtOH_pure) # M_solute missing
    # Edge Cases
    @test cconvert(0.0, "v/v" => "Å^-3"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) == 0.0u"Å^-3"

end

@testitem "NumberDensity Conversions" setup=[CConvert] begin

    N_test = 0.005u"Å^-3" # Approx 8.3 M
    C_equiv = uconvert(u"mol/L", N_test / Na) # ~8.303 M
    # Need rho for this molarity, roughly rho_sol_10M=0.90 is closest defined
    rho_sol_Ntest = rho_sol_10M

    # --- Identity ---
    @test cconvert(N_test, "Å^-3" => "Å^-3") == N_test

    # --- NumberDensity to Molarity ---
    @test cconvert(N_test, "Å^-3" => "M") ≈ C_equiv  rtol=1e-6
    @test cconvert(N_test, "Å^-3" => "M") ≈ 8.303u"mol/L"  rtol=1e-4
    # Missing Kwargs (none needed)
    @test cconvert(N_test, "Å^-3" => "M") isa Quantity
    # Edge Cases
    @test cconvert(0.0u"Å^-3", "Å^-3" => "M") == 0.0u"mol/L"

    # --- NumberDensity to Molality --- (Intermediate Molarity)
    # Ethanol case (N_test -> ~16.8 m)
    @test cconvert(N_test, "Å^-3" => "mol/kg"; M_solute=M_EtOH, rho_solution=rho_sol_Ntest) ≈ 17.95u"mol/kg" rtol=1e-3
    # Missing Kwargs (cascades)
    @test_throws UndefKeywordError cconvert(N_test, "Å^-3" => "mol/kg"; M_solute=M_EtOH) # Missing rho_solution
    @test_throws UndefKeywordError cconvert(N_test, "Å^-3" => "mol/kg"; rho_solution=rho_sol_Ntest) # Missing M_solute
    # Edge Cases
    @test cconvert(0.0u"Å^-3", "Å^-3" => "mol/kg"; M_solute=M_EtOH, rho_solution=rho_H2O_pure) == 0.0u"mol/kg"

    # --- NumberDensity to MoleFraction --- (Intermediate Molarity)
    # Ethanol case (N_test -> χ ~ 0.268)
    @test cconvert(N_test, "Å^-3" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_Ntest) ≈ 0.2443 rtol=1e-3
     # Missing Kwargs (cascades)
    @test_throws UndefKeywordError cconvert(N_test, "Å^-3" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) # Missing rho_solution
    # Edge Cases
    @test cconvert(0.0u"Å^-3", "Å^-3" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_H2O_pure) == 0.0

    # --- NumberDensity to MassFraction --- (Intermediate Molarity)
    # Ethanol case (N_test -> mf ~ 0.425)
    @test cconvert(N_test, "Å^-3" => "w/w"; M_solute=M_EtOH, rho_solution=rho_sol_Ntest) ≈ 0.452649 rtol=1e-3
    # Output Percentage
    @test cconvert(N_test, "Å^-3" => "%m/m"; M_solute=M_EtOH, rho_solution=rho_sol_Ntest) ≈ 45.26 rtol=1e-3
    # Missing Kwargs (cascades)
    @test_throws UndefKeywordError cconvert(N_test, "Å^-3" => "w/w"; M_solute=M_EtOH) # Missing rho_solution
     # Edge Cases
    @test cconvert(0.0u"Å^-3", "Å^-3" => "w/w"; M_solute=M_EtOH, rho_solution=rho_H2O_pure) == 0.0

    # --- NumberDensity to VolumeFraction --- (Intermediate Molarity)
    # Ethanol case (N_test -> vf ~ 0.485)
    @test cconvert(N_test, "Å^-3" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 0.4872 rtol=1e-3
    # Output Percentage
    @test cconvert(N_test, "Å^-3" => "%v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) ≈ 48.72 rtol=1e-3
    # Missing Kwargs (cascades)
    @test_throws UndefKeywordError cconvert(N_test, "Å^-3" => "v/v"; M_solute=M_EtOH) # Missing rho_solute
    # Edge Cases
    @test cconvert(0.0u"Å^-3", "Å^-3" => "v/v"; M_solute=M_EtOH, rho_solute=rho_EtOH_pure) == 0.0

end

@testitem "Wrapper Functionality" setup=[CConvert] begin

    # Test % handling
    @test cconvert(50.0, "%m/m" => "w/w") == 0.5
    @test cconvert(0.5, "w/w" => "%m/m") == 50.0
    @test cconvert(50.0, "%m/m" => "MassFraction") == 0.5 # Target type name
    @test cconvert(0.5, "MassFraction" => "%m/m") == 50.0 # Source type name

    # Error fraction input > 1
    @test_throws ArgumentError cconvert(1.1, "w/w" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O)
    @test cconvert(50.0, " %m/m " => " Molarity "; M_solute=M_EtOH, rho_solution=rho_sol_50ww) ≈ 9.920u"mol/L" rtol=1e-3

    # Test unknown units
    @test_throws ArgumentError cconvert(1.0, "bad_unit" => "Molarity")
    @test_throws ArgumentError cconvert(1.0, "Molarity" => "bad_unit")
    @test_throws ArgumentError cconvert(1.0u"mol/kg", "Molarity" => "Molarity")

    # Test missing kwargs via wrapper
    @test_throws UndefKeywordError cconvert(0.5, "w/w" => "Molarity"; M_solute=M_EtOH) # Missing rho_solution

end

@testitem "Round Trips" setup=[CConvert] begin

    # Use χ=0.2 point
    chi1 = 0.2
    M1 = cconvert(chi1, "chi" => "M"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_chi02)
    m1 = cconvert(chi1, "chi" => "mol/kg"; M_solvent=M_H2O)
    mf1_perc = cconvert(chi1, "chi" => "%m/m"; M_solute=M_EtOH, M_solvent=M_H2O)
    vf1_perc = cconvert(chi1, "chi" => "%v/v"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_chi02)
    N1 = cconvert(chi1, "chi" => "Å^-3"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_chi02)

    # Test round trips back to mole fraction (expect approx equality)
    @test cconvert(M1, "M" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_chi02) ≈ chi1 rtol=1e-3
    @test cconvert(m1, "mol/kg" => "chi"; M_solvent=M_H2O) ≈ chi1 rtol=1e-3
    @test cconvert(mf1_perc, "%m/m" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O) ≈ chi1 rtol=1e-3
    @test cconvert(vf1_perc, "%v/v" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solute=rho_EtOH_pure, rho_solution=rho_sol_chi02) ≈ chi1 rtol=1e-3
    @test cconvert(N1, "Å^-3" => "chi"; M_solute=M_EtOH, M_solvent=M_H2O, rho_solution=rho_sol_chi02) ≈ chi1 rtol=1e-3

end
