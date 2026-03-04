#
# Pre-optimization: adjust molecule positions to satisfy geometric constraints
# before running the main packing optimization with distance penalties.
#
# The algorithm (following the original Fortran Packmol):
# 1. Run a short constraint-only optimization (SPGBox, few iterations)
# 2. Identify "bad" molecules: those with non-zero constraint penalty or
#    overlapping with fixed atoms
# 3. Randomly re-place a fraction (`movefrac`) of the bad molecules
#    within the per-type CM bounds (if provided)
# 4. Repeat up to `nloop` times until all constraints are satisfied.
#
function adjust_constraints!(
    packmol_system::PackmolSystem{D,T},
    free_mol_indices::Vector{Int},
    RNG;
    nloop::Int=200,
    movefrac::T=T(0.05),
    opt_nit::Int=20,
    precision::T=T(1e-10),
    iprint::Int=10,
    domovebad::Bool=true,
    cm_lo_type::Union{Nothing,Vector{SVector{D,T}}} = nothing,
    cm_hi_type::Union{Nothing,Vector{SVector{D,T}}} = nothing,
    buffers::Union{Nothing,MemoryBuffers} = nothing,
) where {D,T}
    nfree = length(free_mol_indices)
    nfree == 0 && return packmol_system

    # Check if there are any constraints to satisfy
    has_constraints = false
    for st in packmol_system.structure_types
        if !isempty(st.constraints)
            has_constraints = true
            break
        end
    end
    has_constraints || return packmol_system

    natoms = length(packmol_system.atoms)
    atom_positions = isnothing(buffers) ? Vector{SVector{D,T}}(undef, natoms) : buffers.atom_positions
    fg_output = InteratomicDistanceFG{D,T}(packmol_system)

    # Build CellListMap system for fixed atom positions (for overlap checks).
    # Only when domovebad=true: at that point molecules are already within their
    # constraint regions. When domovebad=false (initial pass), molecules are still
    # at sidemax (±1000 Å) coordinates; updating fixed_sys with those positions
    # would force CellListMap to build a ~10⁹-cell grid (gigabytes of memory).
    fixed_sys = domovebad ? build_fixed_particle_system(packmol_system) : nothing
    tol = packmol_system.tolerance

    # Build molecule → structure type mapping
    mol_structure_type = _build_mol_structure_type(packmol_system)

    println("  Adjusting initial point to fit the constraints")

    for iloop in 1:nloop
        # Set up optimizer variables from current molecule positions
        x = isnothing(buffers) ? Vector{T}(undef, nfree * 2 * D) : buffers.x
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(free_mol_indices)
            x_mol[k] = packmol_system.molecule_positions[imol]
        end
        auxvecs = isnothing(buffers) ? SPGBox.VAux(x, zero(T)) : buffers.vaux

        # Run a short constraint-only optimization
        spgresult = spgbox!(
            (g, x) -> constraint_only_fg!(g, x, fg_output, packmol_system, atom_positions, free_mol_indices),
            x;
            vaux=auxvecs,
            nitmax=opt_nit,
            nfevalmax=10 * opt_nit,
            callback=(result) -> result.f < precision,
        )

        # Update molecule positions from optimizer
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(free_mol_indices)
            packmol_system.molecule_positions[imol] = x_mol[k]
        end

        # Compute atom positions for per-molecule evaluation
        compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)

        # Compute per-molecule badness (constraint violations + overlap with fixed atoms)
        badness = molecule_badness(packmol_system, atom_positions, fixed_sys, tol; buffers)

        # Check if all constraints are satisfied and no overlaps with fixed
        total_badness = sum(badness[imol] for imol in free_mol_indices)
        if total_badness < precision
            @printf("  Constraint adjustment converged at loop %d (f = %.2e)\n", iloop, total_badness)
            return packmol_system
        end

        # Identify bad molecules among free molecules
        bad_mols = Int[]
        bad_vals = T[]
        for imol in free_mol_indices
            if badness[imol] > precision
                push!(bad_mols, imol)
                push!(bad_vals, badness[imol])
            end
        end
        nbad = length(bad_mols)

        if nbad == 0
            @printf("  All molecules satisfy constraints at loop %d\n", iloop)
            return packmol_system
        end

        # Sort bad molecules from worst to best
        perm = sortperm(bad_vals; rev=true)
        bad_mols = bad_mols[perm]

        # Move a fraction of the bad molecules to new random positions
        nmove = domovebad ? max(1, min(nbad, round(Int, movefrac * nfree))) : 0
        if iloop % iprint == 0 || iloop == 1
            @printf("  Loop %4d: f = %10.4e  bad molecules: %d/%d  moving: %d\n",
                iloop, spgresult.f, nbad, nfree, nmove)
        end
        for i in 1:nmove
            ist = mol_structure_type[bad_mols[i]]
            st = packmol_system.structure_types[ist]
            if !isnothing(cm_lo_type) && !isnothing(cm_hi_type)
                lo = cm_lo_type[ist]
                hi = cm_hi_type[ist]
                has_valid_bounds = all(lo .< hi)
                randomize_molecule!(packmol_system, bad_mols[i], st, RNG;
                    cm_lo = has_valid_bounds ? lo : nothing,
                    cm_hi = has_valid_bounds ? hi : nothing,
                )
            else
                randomize_molecule!(packmol_system, bad_mols[i], st, RNG)
            end
        end
    end

    # Final evaluation
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)
    badness = molecule_badness(packmol_system, atom_positions, fixed_sys, tol; buffers)
    total_badness = sum(badness[imol] for imol in free_mol_indices)
    @printf("  WARNING: constraint adjustment did not fully converge after %d loops (f = %.2e)\n", nloop, total_badness)
    return packmol_system
end
