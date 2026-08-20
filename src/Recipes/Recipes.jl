
using Unitful
using PDBTools: read_pdb, write_pdb, Atom, mass, charge, maxmin, eachresidue

export density_pure_solvent, density_pure_cossolvent
export write_packmol_input
export SolutionBoxUSC
export SolutionBoxUS
export SolutionBoxUWI

const RecipesDirectory = @__DIR__

# Recipes are higher-level, parameter-driven system setups (target densities/
# concentrations instead of explicit molecule counts) built on top of the
# packing engine: solvation boxes today (SolutionBoxUS/USC/UWI), with
# membranes, vesicles, nanotubes, and special (octahedral/dodecahedral) box
# shapes planned. Each recipe implements `write_packmol_input` (and,
# generically below, `packmol`).
abstract type Recipe end

#
# Shared machinery for the `packmol(system::Recipe; ...)` methods (one per
# concrete recipe type, defined alongside each type's `write_packmol_input`
# method): both build a `PackmolSystem` directly in memory — reusing the same
# box-sizing/molecule-count computation as `write_packmol_input` — and hand it
# to `packmol(::PackmolSystem)`, without ever going through a `.inp` file.
#

# The fixed solute, centered at the origin: equivalent to the input file's
# `number 1` / `center` / `fixed 0. 0. 0. 0. 0. 0.` lines (see
# `_apply_fixed_center!` in StructureType.jl for why `center=:geometric`,
# not the `structure_type` default of `:none`, is required to match that
# text-file behavior: `:none` leaves the raw PDB coordinates un-recentered).
_fixed_solute_structure_type(pdbfile::String; tolerance::Real=2.0) = structure_type(
    pdbfile; number=1, tolerance, fixed=(zeros(3), zeros(3)), center=:geometric,
)

# An orthorhombic unit cell of side lengths `2l`, centered at the origin —
# equivalent to the input file's `pbc -l1 -l2 -l3 l1 l2 l3` line.
_recipe_unitcell(l::AbstractVector{<:Quantity}) = (
    unitcell = Matrix{Float64}(Diagonal(2 .* ustrip.(u"Å", l))),
    unitcell_center = zero(SVector{3,Float64}),
)

include("./DensityTable.jl")

#include("concentration.jl")
include("./concentration_units.jl")

# System types
include("./SolutionBoxUS.jl")
include("./SolutionBoxUSC.jl")
include("./SolutionBoxUWI.jl")
