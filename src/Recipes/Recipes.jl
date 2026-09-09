
using Unitful
using PDBTools: read_pdb, write_pdb, Atom, mass, charge, maxmin, eachresidue

export density_pure_solvent, density_pure_cossolvent
export write_packmol_input
export SolutionBoxUSC
export SolutionBoxUS
export SolutionBoxUWI
export Membrane

const RecipesDirectory = @__DIR__

# Recipes are higher-level, parameter-driven system setups (target densities/
# concentrations instead of explicit molecule counts) built on top of the
# packing engine: solvation boxes today (SolutionBoxUS/USC/UWI), each with a
# `pbc` keyword (`:cubic`/`:orthorhombic`/`:dodecahedral`/`:octahedral`)
# selecting the periodic cell shape; membranes, vesicles, and nanotubes are
# planned as future recipes. Each recipe implements `write_packmol_input`
# (and, generically below, `packmol`).
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

# A recipe's periodic unit cell (built by `set_unitcell`, in Å), centered at
# the origin (where the solute is fixed) — equivalent to the input file's
# `unitcell a b c α β γ` line.
_recipe_unitcell(unitcell::AbstractMatrix{Float64}) = (
    unitcell = Matrix{Float64}(unitcell),
    unitcell_center = zero(SVector{3,Float64}),
)

# One-line human-readable description of a recipe's unit cell, for the
# printed/written summary (`a=... Å, b=... Å, c=... Å, α=...°, β=...°, γ=...°`).
function _unitcell_description(unitcell::AbstractMatrix{Float64})
    a, b, c, α, β, γ = _unitcell_abc_angles(unitcell)
    return "a=$a Å, b=$b Å, c=$c Å, α=$(α)°, β=$(β)°, γ=$(γ)°"
end

include("./DensityTable.jl")

#include("concentration.jl")
include("./concentration_units.jl")

# System types
include("./SolutionBoxUS.jl")
include("./SolutionBoxUSC.jl")
include("./SolutionBoxUWI.jl")
include("./Membrane.jl")
