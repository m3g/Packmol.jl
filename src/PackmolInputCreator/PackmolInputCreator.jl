
using PDBTools: read_pdb, write_pdb, Atom, mass, charge, maxmin, eachresidue

export convert_concentration, convert_density_table!
export density_pure_solvent, density_pure_cossolvent
export write_packmol_input
export SolutionBoxUSC
export SolutionBoxUS
export SolutionBoxUWI

abstract type SolutionBox end
density_pure_solvent(system::SolutionBox) = system.density_table[begin, 2]
density_pure_cossolvent(system::SolutionBox) = system.density_table[end, 2]

const PackmolInputCreatorDirectory = @__DIR__

function unit_name(u::String)
    u == "mol/L" && return "molarity"
    u == "x" && return "molar fraction"
    u == "vv" && return "volume fraction"
    u == "mm" && return "mass fraction"
end

# conversion factor from mL/mol to Å^3/molecule
const CMV = 1e24 / 6.02e23

# conversion factor from mol/L to molecules/Å^3
const CMC = 6.02e23 / 1e27

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
include("./SolutionBoxUS.jl")
include("./SolutionBoxUSC.jl")
include("./SolutionBoxUWI.jl")

# Tests
include("./test/runtests.jl")


