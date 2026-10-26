module Packmol

using TestItems: @testitem, @testsnippet
using StaticArrays: SVector, SMatrix, @SMatrix, MMatrix
using LinearAlgebra: norm, eigen, dot, det, Diagonal, svd, cross
using Statistics: mean
using Base: @kwdef
using Base.Threads: @spawn
using ChunkSplitters: chunks, index_chunks
using PDBToolsBase: Atom, center_of_mass, read_pdb, write_pdb
using SPGBox
using Printf: @printf, @sprintf
import CellListMap
# `get_atoms` is a shared function name (defined in PDBToolsBase): imported,
# not defined here, so that it and other packages in the M3G ecosystem
# (MolSimToolkit.jl, ...) can each add their own method to the very same
# generic function instead of each declaring an unrelated one that happens
# to share a name.
import PDBToolsBase: get_atoms

const src_dir = @__DIR__

# API: exported functions
export read_packmol_input
export write_output
export packmol
export PackmolSystem
export StructureType
export structure_type
export get_atoms

# Terminal input/output, progress bar and task synchronization (compatible with `juliac --trim`)
include("./io/terminal.jl")
include("./io/progress_bar.jl")
include("./threads.jl")

# Constraints
include("./constraints/constraints_base.jl")
include("./constraints/boxes.jl")
include("./constraints/spheres.jl")
include("./constraints/planes.jl")
include("./constraints/cylinders.jl")
include("./constraints/ellipsoids.jl")
include("./constraints/gaussians.jl")
include("./constraints/waves.jl")
include("./constraints/exponentials.jl")
include("./constraints/constraint_types.jl")

# Data structures
include("./data_structures/atoms_and_molecules.jl")
include("./data_structures/StructureType.jl")
include("./data_structures/PackmolSystem.jl")
include("./data_structures/periodic_cells.jl")
include("./data_structures/dodecahedron.jl")
include("./data_structures/octahedron.jl")
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

# Input and output
include("./io/pick_file.jl")
include("./io/write_output.jl")

# Packing functions
include("./packmol_main.jl")

# Runner for the legacy packmol (methods in the PackmolJLLExt extension)
include("./packmol_runner.jl")
include("./CLI.jl")

# Packmol input file creator
include("./Recipes/Recipes.jl")

end
