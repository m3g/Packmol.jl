module Packmol

using TestItems: @testitem
using StaticArrays: SVector, SMatrix, @SMatrix
using LinearAlgebra: norm
using Statistics: mean
using Base: @kwdef
using Base.Threads: @spawn
using PDBTools: readPDB, coor 

const src_dir = @__DIR__

# API: exported functions
export pack_monoatomic!

# Constraints
include("./constraints/constraints_base.jl")
include("./constraints/boxes.jl")
include("./constraints/spheres.jl")

# Data structures
include("./data_structures/StructureType.jl")
include("./data_structures/PackmolSystem.jl")

# Read input files


# Rigid body transformations 
include("./rigid_body/rigid_body.jl")
include("./rigid_body/chain_rule.jl")

# Monatomic packing
include("./mono_atomic.jl")

# Runner for the legacy packmol 
include("./packmol_runner.jl")

end
