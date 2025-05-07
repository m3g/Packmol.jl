
using Unitful
using PDBTools: read_pdb, write_pdb, Atom, mass, charge, maxmin, eachresidue

export convert_concentration, convert_density_table!
export density_pure_solvent, density_pure_cossolvent
export write_packmol_input
export SolutionBoxUSC
export SolutionBoxUS
export SolutionBoxUWI

abstract type SolutionBox end
density_pure_solvent(system::SolutionBox) = system.density_table[begin, 2]
# Will be the pure cosolvent if the density table is complete
density_highest_cosolvent_concentration(system::SolutionBox) = system.density_table[end, 2] 

const PackmolInputCreatorDirectory = @__DIR__

#include("concentration.jl")
include("./concentration_units.jl")


# System types
include("./SolutionBoxUS.jl")
#include("./SolutionBoxUSC.jl")
#include("./SolutionBoxUWI.jl")

# Tests
#include("./test/runtests.jl")


