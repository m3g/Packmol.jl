module Packmol

using TestItems: @testitem, @testsnippet
using StaticArrays: SVector, SMatrix, @SMatrix, MMatrix
using LinearAlgebra: norm, eigen, dot, det, Diagonal
using Statistics: mean
using Base: @kwdef
using Base.Threads: @spawn
using ChunkSplitters: chunks, index_chunks
using PDBTools: Atom, read_pdb, write_pdb, coor, center_of_mass
using ProgressMeter: Progress, next!, finish!
using SPGBox
using Printf: @printf, @sprintf
import CellListMap
import Scratch
import Libdl

const src_dir = @__DIR__

# API: exported functions
export read_packmol_input
export write_output
export packmol
export PackmolSystem
export StructureType
export structure_type

# Constraints
include("./constraints/constraints_base.jl")
include("./constraints/boxes.jl")
include("./constraints/spheres.jl")
include("./constraints/planes.jl")
include("./constraints/cylinders.jl")
include("./constraints/ellipsoids.jl")
include("./constraints/constraint_types.jl")

# Data structures
include("./data_structures/atoms_and_molecules.jl")
include("./data_structures/StructureType.jl")
include("./data_structures/PackmolSystem.jl")
include("./data_structures/FixedParticleSystem.jl")

# Random number generation
include("./initial_approximation/random.jl")

# Rigid body transformations
include("./rigid_body/rigid_body.jl")
include("./rigid_body/chain_rule.jl")

# Function and gradient of the distance between atoms
include("./function_and_gradient/interatomic_distance_fg.jl")

# Initial approximation: placement, constraint pre-optimization
include("./initial_approximation/initial_approximation.jl")

# Output
include("./io/write_output.jl")

# GENCAN FFI wrapper (experimental, opt-in alternative optimizer)
include("./gencan/build.jl")
include("./gencan/ffi.jl")
include("./gencan/gencan.jl")

# Packing functions
include("./packmol_main.jl")

# Runner for the legacy packmol 
include("./packmol_runner.jl")
@static if VERSION >= v"1.12" 
    include("./CLI.jl")
end

# Packmol input file creator
include("./Recipes/Recipes.jl")

function __init__()
    # A @cfunction pointer is a JIT address valid only for the current
    # process; it must be computed here, not as a top-level `const`, or it
    # would be a stale/dangling pointer after loading a precompiled image.
    GENCAN_FG_CALLBACK[] = @cfunction(
        gencan_fg_callback, Cint, (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble})
    )
    return nothing
end

end
