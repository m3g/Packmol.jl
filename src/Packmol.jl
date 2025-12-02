module Packmol

using TestItems: @testitem, @testsnippet
using StaticArrays: SVector, SMatrix, @SMatrix, MMatrix
using LinearAlgebra: norm, eigen
using Statistics: mean
using Base: @kwdef
using Base.Threads: @spawn
using PDBTools: Atom, readPDB, coor 
import CellListMap

const src_dir = @__DIR__

# API: exported functions
export pack_monoatomic!

# Constraints
include("./constraints/constraints_base.jl")
include("./constraints/boxes.jl")
include("./constraints/spheres.jl")

# Data structures
include("./data_structures/atoms_and_molecules.jl")
include("./data_structures/StructureType.jl")
include("./data_structures/PackmolSystem.jl")

# Read input files


# Rigid body transformations 
include("./rigid_body/rigid_body.jl")
include("./rigid_body/chain_rule.jl")

# Monatomic packing
include("./mono_atomic.jl")

# Function and gradient of the distance between atoms
include("./interatomic_distance_fg.jl")

# Runner for the legacy packmol 
include("./packmol_runner.jl")
@static if VERSION >= v"1.12" 
    include("./CLI.jl")
end

end
