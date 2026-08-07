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

    tol = packmol_system.tolerance
    has_pbc = !isnothing(packmol_system.unitcell)
    max_natoms = maximum(st.natoms for st in packmol_system.structure_types if !st.fixed.fixed; init=0)
    ntasks = Threads.nthreads()

    max_extent = zero(T)
    for st in packmol_system.structure_types
        st.fixed.fixed && continue
        for r in st.reference_coordinates
            max_extent = max(max_extent, norm(r))
        end
    end

    # Build the fixed-atom overlap-check system once, and deep-copy it once per
    # thread — not once per structure type, and not once per trial molecule.
    #
    # Two things used to make this step slow and memory-hungry on systems
    # with a large fixed structure (e.g. a membrane/capsid) and a free-molecule
    # placement region much bigger than that structure's own footprint:
    #
    # 1. With no explicit unitcell, CellListMap treats the system as
    #    non-periodic: every update!/pairwise! call — one per trial molecule,
    #    up to `nfree * max_guess_try` times — recomputes
    #    `limits(xpositions, ypositions)`, an O(n_fixed) scan over *all* fixed
    #    atoms, and rebuilds the box (and the fixed-atom cell list with it)
    #    whenever a trial point falls outside it. Giving the system an
    #    explicit box up front makes CellListMap treat it as periodic instead
    #    and skip all of that.
    # 2. That box must not simply be inflated to cover the whole free-molecule
    #    placement region: cell count scales with box volume / cutoff^D, and
    #    the placement region can be far larger than the fixed structure
    #    itself, which would blow up the cell list even though most of that
    #    volume has no fixed atoms in it at all. Instead the box is sized to
    #    just the fixed atoms' own extent (+ margin), and overlaps_fixed
    #    rejects trial molecules whose bounding box can't possibly reach the
    #    fixed atoms *before* calling into CellListMap at all — see there.
    fixed_atom_positions = _collect_fixed_positions(packmol_system)
    fixed_sys, fixed_lo, fixed_hi = if isempty(fixed_atom_positions)
        nothing, zero(SVector{D,T}), zero(SVector{D,T})
    else
        lo, hi = compute_bounding_box(fixed_atom_positions)
        margin = max_extent + tol
        fixed_unitcell = T(1.2) * ((hi .+ margin) - (lo .- margin))
        volume = prod(fixed_unitcell)
        cutoff = _capped_cutoff(volume, tol, length(fixed_atom_positions), D)
        sys = build_fixed_particle_system(packmol_system; nmols=1, parallel=false, unitcell=fixed_unitcell, cutoff)
        sys, lo, hi
    end
    task_fixed_syss = isnothing(fixed_sys) ? nothing : [deepcopy(fixed_sys) for _ in 1:ntasks]

    println("  Setting initial trial coordinates (best of $max_guess_try per molecule)...")

    # Structure types are processed sequentially, each under its own @sync, so
    # that reusing task_fixed_syss[ichunk] across structure types is race-free:
    # all of type A's tasks finish before type B's tasks — which touch the
    # same objects — are spawned.
    imol_offset = 0
    for (ist, st) in enumerate(packmol_system.structure_types)
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

        @sync for (ichunk, irange) in enumerate(chunks(1:st.number_of_molecules; n=ntasks))
            task_seed = rand(RNG, UInt64)
            task_fixed_sys = isnothing(task_fixed_syss) ? nothing : task_fixed_syss[ichunk]
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
                           overlaps_fixed(mp, st.reference_coordinates, task_fixed_sys, mol_positions, fixed_lo, fixed_hi, tol)
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
