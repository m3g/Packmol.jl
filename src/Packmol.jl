module Packmol

using TestItems: @testitem, @testsnippet
using StaticArrays: SVector, SMatrix, @SMatrix, MMatrix
using LinearAlgebra: norm, eigen, dot, det, Diagonal
using Statistics: mean
using Base: @kwdef
using Base.Threads: @spawn
using ChunkSplitters: chunks
using PDBTools: Atom, readPDB, writePDB, coor, center_of_mass
import CellListMap
using SPGBox
using Printf: @printf, @sprintf

const src_dir = @__DIR__

# API: exported functions
export read_packmol_input
export write_output
export packmol

# Constraints
include("./constraints/constraints_base.jl")
include("./constraints/boxes.jl")
include("./constraints/spheres.jl")

# Data structures
include("./data_structures/atoms_and_molecules.jl")
include("./data_structures/StructureType.jl")
include("./data_structures/PackmolSystem.jl")

# Random number generation
include("./random.jl")

# Rigid body transformations
include("./rigid_body/rigid_body.jl")
include("./rigid_body/chain_rule.jl")

# Monatomic packing
include("./mono_atomic.jl")

# Function and gradient of the distance between atoms
include("./interatomic_distance_fg.jl")

# Initial approximation: placement, constraint pre-optimization
include("./initial_approximation.jl")

# Runner for the legacy packmol 
include("./packmol_runner.jl")
@static if VERSION >= v"1.12" 
    include("./CLI.jl")
end

# Packmol input file creator
include("./PackmolInputCreator/PackmolInputCreator.jl")

end
