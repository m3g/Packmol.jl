#
# Re-initialize molecule positions using best-of-N random placement within
# per-type CM bounds (following Fortran Packmol initial.f90).
#
# After constraint fitting (adjust_constraints!), the molecule CMs define
# per-type bounding regions. This function generates `max_guess_try` random
# placements for each molecule within those bounds, evaluates the constraint
# penalty for each, and keeps the best placement.
#
function reinitialize_with_bounds!(
    packmol_system::PackmolSystem{D,T},
    free_mol_indices::Vector{Int},
    RNG;
    max_guess_try::Int=20,
    precision::T=T(1e-2),
    cm_min::Union{Nothing,Vector{SVector{D,T}}} = nothing,
    cm_max::Union{Nothing,Vector{SVector{D,T}}} = nothing,
) where {D,T}
    nfree = length(free_mol_indices)
    nfree == 0 && return packmol_system

    # Use provided CM bounds, or compute from current molecule positions
    if isnothing(cm_min) || isnothing(cm_max)
        cm_min, cm_max = compute_cm_bounds(packmol_system)
    end

    # Build fixed particle system for overlap checks
    fixed_sys = build_fixed_particle_system(packmol_system)
    tol = packmol_system.tolerance

    max_natoms = maximum(st.natoms for st in packmol_system.structure_types if !st.fixed.fixed; init=0)

    has_pbc = !isnothing(packmol_system.unitcell)

    println("  Setting initial trial coordinates (best of $max_guess_try per molecule)...")

    imol_offset = 0
    @sync for (ist, st) in enumerate(packmol_system.structure_types)
        st_offset = imol_offset
        if st.fixed.fixed
            imol_offset += st.number_of_molecules
            continue
        end

        # Per-type CM bounds from constraint fitting
        lo = cm_min[ist]
        hi = cm_max[ist]
        extent = hi - lo

        # Fallback: if bounds are degenerate (single molecule or all at same CM),
        # use PBC unit cell or sidemax box
        if any(extent .≤ zero(T))
            if has_pbc
                uc = packmol_system.unitcell
                center = packmol_system.unitcell_center
            else
                sidemax = T(DEFAULT_SIDEMAX)
                lo = SVector{D,T}(ntuple(_ -> -sidemax, D))
                hi = SVector{D,T}(ntuple(_ -> sidemax, D))
                extent = hi - lo
            end
        end

        for (_, irange) in enumerate(chunks(1:st.number_of_molecules; n=Threads.nthreads()))
            task_seed = rand(RNG, UInt64)
            # Each task gets its own copy of the fixed system so that concurrent
            # update!/pairwise! calls don't race. The fixed-atom (y) cell list is
            # already built inside fixed_sys and is deep-copied once per task.
            task_fixed_sys = isnothing(fixed_sys) ? nothing : deepcopy(fixed_sys)
            @spawn begin
                task_rng = typeof(RNG)(task_seed)
                mol_positions = Vector{SVector{D,T}}(undef, max_natoms)
                for i in irange
                    imol = st_offset + i
                    best_mp = packmol_system.molecule_positions[imol]
                    # Use the actual step-2 penalty as baseline: only replace if
                    # a random trial is strictly better. If the molecule already
                    # satisfies its constraints, skip the trial loop entirely —
                    # don't undo the constraint optimization's work.
                    best_fx = constraint_penalty_sum(best_mp, st)
                    if best_fx < precision
                        continue
                    end

                    for itry in 1:max_guess_try
                        # Random CM within per-type bounds
                        if has_pbc && any(extent .≤ zero(T))
                            # PBC fallback: random within unit cell
                            frac = SVector{D,T}(ntuple(_ -> rand(task_rng, T) - T(0.5), D))
                            cm = SVector{D,T}(uc * frac) + center
                        else
                            cm = lo + SVector{D,T}(ntuple(d -> rand(task_rng, T) * extent[d], D))
                        end
                        angles = SVector{D,T}(ntuple(_ -> T(2π) * rand(task_rng, T), D))
                        mp = MoleculePosition(cm, angles)

                        # Check overlap with fixed atoms first (skip if overlapping)
                        if packmol_system.avoid_overlap &&
                           overlaps_fixed(mp, st.reference_coordinates, task_fixed_sys, mol_positions)
                            continue
                        end

                        # Evaluate constraint penalty
                        fx = constraint_penalty_sum(mp, st)
                        if fx < best_fx
                            best_fx = fx
                            best_mp = mp
                        end

                        # Stop early if constraints satisfied
                        if best_fx < precision
                            break
                        end
                    end
                    packmol_system.molecule_positions[imol] = best_mp
                end
            end
        end
        imol_offset += st.number_of_molecules
    end

    return packmol_system
end
