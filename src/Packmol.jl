#
# The module: its dependencies, the exported names, and the included files, in an order
# in which every type is defined before it is used in a method signature (see README.md).
#
module Packmol

using TestItems: @testitem, @testsnippet
using StaticArrays: SVector, SMatrix, @SMatrix, MMatrix, MVector
using LinearAlgebra: norm, eigen, eigvals, dot, det, Diagonal, svd, cross, Symmetric
using Statistics: mean
using Base: @kwdef
using Base.Threads: @spawn
using ChunkSplitters: chunks, index_chunks
using PDBToolsBase: Atom, Structure, center_of_mass, read_pdb, write_pdb
using PDBToolsBase: read_xyz, write_xyz, read_mmcif
using SPGBox
using Printf: @printf, @sprintf
import CellListMap
using CellListMap: ParticleSystem, NeighborPair, pairwise!
using FunctionWrappers: FunctionWrapper
import Random
# `get_atoms` is a shared function name (defined in PDBToolsBase): imported,
# not defined here, so that it and other packages in the M3G ecosystem
# (MolSimToolkit.jl, ...) can each add their own method to the very same
# generic function instead of each declaring an unrelated one that happens
# to share a name.
import PDBToolsBase
import PDBToolsBase: get_atoms

const src_dir = @__DIR__

# Set when the `packmol` executable is built (juliac/build.jl): the precompilation
# workload is skipped, and the constraints available only in Julia code (CustomConstraint)
# are left out of the constraint union (see constraints/constraint_types.jl), since
# `juliac --trim` does not compile them.
const _executable_build = get(ENV, "PACKMOL_EXECUTABLE_BUILD", "false") == "true"

# test_file(name): path of the test file `name` (e.g. "structure_files/water.pdb"). The
# test files are a lazy artifact (see Artifacts.toml), downloaded on the first call. The
# method is defined in the TestFilesExt extension, available after `using LazyArtifacts`.
function test_file end

# API: exported functions
export read_packmol_input
export write_output
export packmol
export guess_volume, free_volume
export voxel_place
export PackmolSystem
export StructureType
export structure_type
export get_atoms
export shape_descriptors, select_conformations

# Basic utilities: terminal output, progress bar and task synchronization (compatible
# with `juliac --trim`), and the file dialog
include("./utils/terminal.jl")
include("./utils/progress_bar.jl")
include("./utils/threads.jl")
include("./utils/pick_file.jl")

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
include("./constraints/custom.jl")
include("./constraints/surfaces.jl")
include("./constraints/constraint_types.jl")
include("./constraints/rotations.jl")

# Data structures
include("./data_structures/atoms_and_molecules.jl")
include("./data_structures/StructureType.jl")
include("./data_structures/PackmolSystem.jl")
include("./data_structures/system_options.jl")
include("./data_structures/FixedParticleSystem.jl")

# Periodic cells
include("./periodic_cells/periodic_cells.jl")
include("./periodic_cells/dodecahedron.jl")
include("./periodic_cells/octahedron.jl")
include("./periodic_cells/wrapping.jl")

# Input files
include("./io/read_input.jl")
include("./io/read_structure_block.jl")

# Rigid body transformations
include("./rigid_body/rigid_body.jl")
include("./rigid_body/chain_rule.jl")
include("./rigid_body/rotation_metric.jl")

# The objective function and its gradient
include("./function_and_gradient/fg_output.jl")
include("./function_and_gradient/pair_fg.jl")
include("./function_and_gradient/atom_positions.jl")
include("./function_and_gradient/constraint_fg.jl")
include("./function_and_gradient/fg.jl")

# Initial placement of the molecules, and movebad
include("./initial_approximation/bounds.jl")
include("./initial_approximation/fixed_molecules.jl")
include("./initial_approximation/constraint_checks.jl")
include("./initial_approximation/buffers.jl")
include("./initial_approximation/movebad.jl")
include("./initial_approximation/initialize_molecules.jl")
include("./initial_approximation/molecule_badness.jl")
include("./initial_approximation/randomize_molecule.jl")
include("./initial_approximation/reinitialize_with_bounds.jl")
include("./initial_approximation/adjust_constraints.jl")
include("./initial_approximation/initial_approximation.jl")

# Output and restart
include("./io/write_output.jl")
include("./io/restart.jl")

# The packing
include("./packing/packmol.jl")
include("./packing/packing_loop.jl")
include("./packing/report.jl")
include("./packing/active_set.jl")

# Tools for setting up systems: volumes, shapes of molecules, and tools for building
# systems by hand (public, not exported)
include("./tools/volume.jl")
include("./tools/shape.jl")
include("./tools/system_building.jl")
include("./tools/voxel_place.jl")

# Runner for the legacy packmol (methods in the PackmolJLLExt extension)
include("./packmol_runner.jl")

# Command line interface
include("./CLI.jl")

# Recipes
include("./Recipes/Recipes.jl")

# Precompilation workload
include("./precompile.jl")

end
