
using Unitful
using PDBTools: read_pdb, write_pdb, Atom, mass, charge, maxmin, eachresidue

export convert_concentration, convert_density_table!
export density_pure_solvent, density_pure_cossolvent
export write_packmol_input
export SolutionBoxUSC
export SolutionBoxUS
export SolutionBoxUWI

const PackmolInputCreatorDirectory = @__DIR__

abstract type SolutionBox end

include("./DensityTable.jl")

#include("concentration.jl")
include("./concentration_units.jl")


# System types
include("./SolutionBoxUS.jl")
include("./SolutionBoxUSC.jl")
include("./SolutionBoxUWI.jl")
