using Unitful
using OrderedCollections: OrderedDict
using Compat: @compat

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

# Round for error checking, to avoid floating point issues
_round(x::Real) = round(x; digits=3)
_round(x::T) where {T<:Quantity} = round(T, x; digits=3)

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

    _round(m_solute) >= _round(m_solution) && throw(ArgumentError("Solute mass ($m_solute) >= solution mass ($m_solution)"))
    m_solvent = m_solution - m_solute
    _round(m_solvent) < 0 * unit(m_solvent) && throw(ArgumentError("Non-positive solvent mass ($m_solvent)"))

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

    _round(m_solute) > _round(m_solution) && throw(ArgumentError("Solute mass ($m_solute) >= solution mass ($m_solution)"))
    m_solvent = m_solution - m_solute
    _round(m_solvent) < 0 * unit(m_solvent) && throw(ArgumentError("Non-positive solvent mass ($m_solvent)"))

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

    _round(m_solution) < 0 * unit(m_solution) && throw(ArgumentError("Non-positive solution mass ($m_solution)"))
    _round(m_solute) > _round(m_solution) && throw(ArgumentError("Solute mass ($m_solute) > solution mass ($m_solution)"))

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

    _round(m_solution) < 0 * unit(m_solution) && throw(ArgumentError("Non-positive solution mass ($m_solution)"))
    V_solution = uconvert(DEFAULT_VOLUME_UNIT, m_solution / rho_s)
    _round(V_solution) < 0 * unit(V_solution) && throw(ArgumentError("Non-positive solution volume ($V_solution)"))

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

    _round(m_solution) < 0 * unit(m_solution) && throw(ArgumentError("Non-positive solution mass ($m_solution)"))

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
    !(0 <= _round(chi_val) <= 1) && throw(ArgumentError("Mole fraction must be 0-1"))
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

    _round(m_solution) < 0 * unit(m_solution) && throw(ArgumentError("Non-positive solution mass ($m_solution)"))
    V_solution = uconvert(DEFAULT_VOLUME_UNIT, m_solution / rho_s)
    _round(V_solution) < 0 * unit(V_solution) && throw(ArgumentError("Non-positive solution volume ($V_solution)"))

    molarity = n_solute / V_solution
    return uconvert(DEFAULT_MOLARITY_UNIT, molarity)
end

function cconvert(
    chi_val::Real, ::Type{MoleFraction}, ::Type{Molality};
    M_solvent, kwargs...
)
    !(0 <= _round(chi_val) < 1) &&  throw(ArgumentError("Mole fraction must be 0 <= χ < 1 for Molality"))
    chi = chi_val * NoUnits
    Msv = _ensure_unit(M_solvent, DEFAULT_MASS_UNIT)

    n_total_basis = 1u"mol"
    n_solute = chi * n_total_basis
    n_solvent = (1*NoUnits - chi) * n_total_basis

    m_solvent = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solvent * Msv)
    _round(m_solvent) <= 0 * unit(m_solvent) && throw(ArgumentError("Non-positive solvent mass ($m_solvent), χ=1?"))

    molality = n_solute / m_solvent
    return uconvert(DEFAULT_MOLALITY_UNIT, molality)
end

function cconvert(
    chi_val::Real, ::Type{MoleFraction}, ::Type{MassFraction};
    M_solute, M_solvent, kwargs...
)
    !(0 <= _round(chi_val) <= 1) && throw(ArgumentError("Mole fraction must be 0-1"))
    chi = chi_val * NoUnits
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)
    Msv = _ensure_unit(M_solvent, DEFAULT_MASS_UNIT)

    n_total_basis = 1u"mol"
    n_solute = chi * n_total_basis
    n_solvent = (1*NoUnits - chi) * n_total_basis

    m_solute = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solute * Ms)
    m_solvent = uconvert(DEFAULT_MASS_QUANTITY_UNIT, n_solvent * Msv)
    m_solution = m_solute + m_solvent

    _round(m_solution) <= 0 * unit(m_solution) && throw(ArgumentError("Non-positive solution mass ($m_solution)"))

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
    !(0 <= _round(mf_val) <= 1) && throw(ArgumentError("Mass fraction must be 0-1"))
    mf = mf_val * NoUnits
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)
    rho_s = _ensure_unit(rho_solution, DEFAULT_DENSITY_UNIT)

    m_solution_basis = 1u"kg" # Basis mass
    m_solute = mf * m_solution_basis
    n_solute = uconvert(u"mol", m_solute / Ms)

    V_solution = uconvert(DEFAULT_VOLUME_UNIT, m_solution_basis / rho_s)
    _round(V_solution) < 0 * unit(V_solution) && throw(ArgumentError("Non-positive solution volume ($V_solution)"))

    molarity = n_solute / V_solution
    return uconvert(DEFAULT_MOLARITY_UNIT, molarity)
end

function cconvert(
    mf_val::Real, ::Type{MassFraction}, ::Type{Molality};
    M_solute, kwargs...
)
    !(0 <= _round(mf_val) < 1) && throw(ArgumentError("Mass fraction must be 0 <= mf < 1 for Molality"))
    mf = mf_val * NoUnits
    Ms = _ensure_unit(M_solute, DEFAULT_MASS_UNIT)

    m_solution_basis = 1u"kg" # Basis mass
    m_solute = mf * m_solution_basis
    m_solvent = (1*NoUnits - mf) * m_solution_basis

    _round(m_solvent) < 0 * unit(m_solvent) && throw(ArgumentError("Non-positive solvent mass ($m_solvent), mf=1?"))
    n_solute = uconvert(u"mol", m_solute / Ms)

    molality = n_solute / m_solvent
    return uconvert(DEFAULT_MOLALITY_UNIT, molality)
end

function cconvert(
    mf_val::Real, ::Type{MassFraction}, ::Type{MoleFraction};
    M_solute, M_solvent, kwargs...
)
    !(0 <= _round(mf_val) <= 1) && throw(ArgumentError("Mass fraction must be 0-1"))
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
    !(0 <= _round(mf_val) <= 1) && throw(ArgumentError("Mass fraction must be 0-1"))
    mf = mf_val * NoUnits
    rho_pure_s = _ensure_unit(rho_solute, DEFAULT_DENSITY_UNIT)
    rho_s = _ensure_unit(rho_solution, DEFAULT_DENSITY_UNIT)

    m_solution_basis = 1u"kg" # Basis mass
    m_solute = mf * m_solution_basis

    V_solute = uconvert(DEFAULT_VOLUME_UNIT, m_solute / rho_pure_s)
    V_solution = uconvert(DEFAULT_VOLUME_UNIT, m_solution_basis / rho_s)
    _round(V_solution) < 0 * unit(V_solution) && throw(ArgumentError("Non-positive solution volume ($V_solution)"))

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
    !(0 <= _round(vf_val) <= 1) && throw(ArgumentError("Volume fraction must be 0-1"))
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
     !(0 <= _round(vf_val) <= 1) && throw(ArgumentError("Volume fraction must be 0-1"))
     vf = vf_val * NoUnits
     rho_pure_s = _ensure_unit(rho_solute, DEFAULT_DENSITY_UNIT)
     rho_s = _ensure_unit(rho_solution, DEFAULT_DENSITY_UNIT)

     V_solution_basis = 1u"L" # Basis volume
     V_solute = vf * V_solution_basis

     m_solute = uconvert(DEFAULT_MASS_QUANTITY_UNIT, V_solute * rho_pure_s)
     m_solution = uconvert(DEFAULT_MASS_QUANTITY_UNIT, V_solution_basis * rho_s)
     _round(m_solution) <= 0 * unit(m_solution) && throw(ArgumentError("Non-positive solution mass ($m_solution)"))

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
@compat public UNIT_STRINGS

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
    using Unitful: ustrip, Na, dimension, Quantity, uconvert

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
