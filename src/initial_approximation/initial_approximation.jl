#
# Functions related to the generation of the initial approximation for the packing problem.
#
# Default sidemax: half-width of the large box for initial placement.
const DEFAULT_SIDEMAX = 1000.0

include("./_clamp_molecules_to_bounds.jl")
include("./compute_bounding_box.jl")
include("./overlaps_fixed.jl")
include("./satisfies_constraints.jl")
include("./constraint_penalty_sum.jl")
include("./compute_cm_bounds.jl")
include("./movebad.jl")
include("./initialize_molecules.jl")
include("./_collect_fixed_positions.jl")
include("./_build_mol_structure_type.jl")
include("./MemoryBuffers.jl")
include("./_fill_free_atoms.jl")
include("./molecule_badness.jl")
include("./randomize_molecule.jl")
include("./reinitialize_with_bounds.jl")
include("./constraint_only_fg.jl")
include("./adjust_constraints.jl")

#
# Set the initial approximation for the packing problem.
#
# Steps:
#   1. Random positions in large sidemax box (initialize_molecules!)
#   2. Constraint-only optimization (no movebad) to push molecules toward
#      feasible regions, then compute CM bounds from the result
#   3. Re-initialize molecules within estimated bounds (best-of-N per molecule)
#   4. Full constraint adjustment with solver and movebad
#   5. Clamp molecules whose CMs are still far outside the constraint bounds
#
function set_initial_approximation!(
    packmol_system::PackmolSystem{D,T},
    free_mol_indices::Vector{Int},
    RNG;
    restart::Bool=false,
    buffers::Union{Nothing,MemoryBuffers{D,T}} = nothing,
) where {D,T}
    # Step 1: Random initial placement
    initialize_molecules!(packmol_system, RNG)

    restart && return packmol_system

    # Step 2: Constraint-only optimization (no movebad) to push molecules
    # toward feasible regions, then compute CM bounds from the result.
    adjust_constraints!(packmol_system, free_mol_indices, RNG; domovebad=false, buffers)
    # Use only constraint-satisfying molecules for bounds: outliers still at
    # sidemax (±1000 Å) would otherwise produce a huge bounding box, causing
    # CellListMap to allocate a ~10⁹-cell grid in subsequent steps.
    cm_min, cm_max = compute_cm_bounds(packmol_system; precision=packmol_system.constraint_precision)

    # Step 3: Re-initialize molecules within the estimated bounds
    # (best-of-N tries per molecule, checking constraints + fixed overlap)
    reinitialize_with_bounds!(packmol_system, free_mol_indices, RNG;
        precision=packmol_system.constraint_precision,
        cm_min=cm_min, cm_max=cm_max,
    )

    # Step 4: Full constraint adjustment with the solver and movebad,
    # randomizing bad molecules within the CM bounds from step 2
    if packmol_system.adjust_constraints_on_init
        adjust_constraints!(packmol_system, free_mol_indices, RNG;
            cm_lo_type=cm_min, cm_hi_type=cm_max, buffers,
        )
    end

    # Step 5: Clamp any molecules whose CMs are still far outside the
    # constraint bounds (can happen when adjustment doesn't converge).
    # This prevents a huge bounding box that would blow up CellListMap memory.
    _clamp_molecules_to_bounds!(packmol_system, free_mol_indices, cm_min, cm_max, RNG)

    return packmol_system
end
