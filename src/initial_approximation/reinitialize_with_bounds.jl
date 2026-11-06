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
    max_guess_try::Int=100,
    precision::T=T(1e-2),
    cm_min::Union{Nothing,Vector{SVector{D,T}}} = nothing,
    cm_max::Union{Nothing,Vector{SVector{D,T}}} = nothing,
    fit::Bool = true,
    max_shift::Union{Nothing,T} = nothing,
    opt_nit::Int = 20,
    # The fits are complete (to the precision of `adjust_constraints!`), not only
    # to `precision`: the molecules left with small violations would be moved
    # randomly by the next adjustment of the constraints
    fit_precision::T = T(1e-10),
) where {D,T}
    nfree = length(free_mol_indices)
    nfree == 0 && return packmol_system

    # Use provided CM bounds, or compute from current molecule positions
    # (`something`, for the bounds to have concrete types)
    cm_lo_type, cm_hi_type = (isnothing(cm_min) || isnothing(cm_max)) ?
        compute_cm_bounds(packmol_system) : (something(cm_min), something(cm_max))

    # Avoid fixed-atom overlaps against the same (loose) distance the main
    # packing loop's own first iteration actually enforces — `radscale *
    # tolerance` (Fortran Packmol's `discale`; see `packmol_system.radscale`'s
    # docstring) — rather than the bare tolerance. Placing free molecules only
    # `tolerance` away from the fixed structure still leaves them violating
    # that first, loosest loop, which just pushes the problem into the
    # (slower) main optimization instead of avoiding it here where it's cheap.
    tol = _global_radscale(packmol_system) * packmol_system.tolerance
    has_pbc, uc, center = _typed_unitcell(packmol_system)
    max_natoms = maximum(st.natoms for st in packmol_system.structure_types if !st.fixed.fixed; init=0)
    ntasks = _nbatches()

    # Build the fixed-atom overlap-check system once per thread (see
    # `_task_overlap_check_systems`) — not once per structure type, and not once per trial molecule.
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
    fixed_sys, fixed_lo, fixed_hi = _build_overlap_check_system(packmol_system, tol)
    task_fixed_syss = _task_overlap_check_systems(fixed_sys, packmol_system, tol, ntasks)

    # (A function barrier: `task_fixed_syss` is `nothing` or a vector of systems, and the
    # closures of the tasks spawned there must have concrete types for `juliac --trim`.)
    _place_molecules!(
        packmol_system, RNG, task_fixed_syss, fixed_lo, fixed_hi, cm_lo_type, cm_hi_type,
        has_pbc, uc, center, tol, precision, max_guess_try, max_natoms, ntasks, nfree,
        fit, something(max_shift, packmol_system.tolerance), opt_nit, fit_precision,
    )
    return packmol_system
end

# Best-of-`max_guess_try` placement of the free molecules (see `reinitialize_with_bounds!`).
function _place_molecules!(
    packmol_system::PackmolSystem{D,T}, RNG, task_fixed_syss, fixed_lo, fixed_hi, cm_min, cm_max,
    has_pbc, uc, center, tol, precision, max_guess_try, max_natoms, ntasks, nfree,
    fit, max_shift, opt_nit, fit_precision,
) where {D,T}
    _println("  Setting initial trial coordinates (best of $max_guess_try per molecule)...")
    progress_meter = ProgressBar(nfree; desc="  Molecules placed: ", barlen=40)

    # Buffers of the constraint fits of the trials, one per task
    natoms = length(packmol_system.atoms)
    task_fg_outputs = [InteratomicDistanceFG{D,T}(packmol_system) for _ in 1:(fit ? ntasks : 0)]
    task_atom_positions = [Vector{SVector{D,T}}(undef, natoms) for _ in 1:(fit ? ntasks : 0)]
    mol_structure_type = _build_mol_structure_type(packmol_system)
    # Bounds of the fit of a trial (see `_fit_bounds`), the same for all
    fit_bounds = _fit_bounds(packmol_system, 1)
    mol_iat_first = _build_mol_iat_first(packmol_system)
    ok_precision = fit ? fit_precision : precision
    task_accepted = zeros(Int, ntasks)
    task_feasible = zeros(Int, ntasks)

    # Structure types are processed sequentially, each under its own @_sync, so
    # that reusing task_fixed_syss[ichunk] across structure types is race-free:
    # all of type A's tasks finish before type B's tasks — which touch the
    # same objects — are spawned.
    wigner_seitz = packmol_system.wigner_seitz_cell
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
        if any(extent .≤ zero(T)) && !has_pbc
            sidemax = T(DEFAULT_SIDEMAX)
            lo = SVector{D,T}(ntuple(_ -> -sidemax, D))
            hi = SVector{D,T}(ntuple(_ -> sidemax, D))
            extent = hi - lo
        end
        # Extended by the radius of the molecule (see `_extended_cm_bounds`)
        if all(extent .> zero(T))
            lo, hi = _extended_cm_bounds(packmol_system, st, lo, hi)
            extent = hi - lo
        end
        # (New names, assigned once, for the tasks below: capturing `lo`/`extent`,
        # which are reassigned above, would box them.)
        cm_lo, cm_extent = lo, extent

        @_sync for (ichunk, irange) in enumerate(chunks(1:st.number_of_molecules; n=ntasks))
            task_seed = rand(RNG, UInt64)
            task_fixed_sys = isnothing(task_fixed_syss) ? nothing : task_fixed_syss[ichunk]
            @spawn begin
                task_rng = typeof(RNG)(task_seed)
                mol_positions = Vector{SVector{D,T}}(undef, max_natoms)
                task_x = Vector{T}(undef, 2 * D)
                task_mol_list = [0]
                for i in irange
                    imol = st_offset + i
                    # Kept only as a last-resort fallback (below) if every
                    # trial overlaps the fixed structure — never trusted as
                    # an automatic "already good enough" bypass. Step 2's own
                    # constraint-only fit pushes every out-of-region molecule
                    # to the *nearest* feasible point, not a random one within
                    # the feasible region: for a thin PBC slab bounded by two
                    # one-sided plane constraints where one side coincides
                    # with the periodic wrap seam (so it's essentially never
                    # violated from a uniform start), that collapses most of
                    # the type onto a single point hugging the other edge —
                    # a real constraint-penalty-satisfying position, but a
                    # badly clustered one, not the uniform-within-bounds
                    # sample this step exists to produce. So every molecule
                    # always gets a fresh random trial below regardless of
                    # how good its inherited position already looks.
                    inherited_mp = packmol_system.molecule_positions[imol]
                    best_mp = inherited_mp
                    # Trials are ranked by (rank, value): rank 0 for those satisfying
                    # the constraints (value: how far the fit moved the center of
                    # mass), 1 for the others (value: the constraint penalty), 2
                    # while there is none (all overlap the fixed structure)
                    best_key = (2, typemax(T))

                    # Two phases: `max_guess_try` trials as drawn, the first satisfying
                    # the constraints accepted (a uniform sample of the region), and only
                    # if none does, `max_guess_try` fitted trials (see below). Fitting
                    # every trial accepted those drawn up to `max_shift` outside the
                    # region, moved onto its surface: 43% more water than uniform within
                    # 3 Å of the outer surface of a vesicle (`outside sphere`).
                    for itry in 1:(fit ? 2 : 1) * max_guess_try
                        # Random CM within per-type bounds
                        if has_pbc && any(cm_extent .≤ zero(T))
                            # PBC fallback: random within unit cell
                            frac = SVector{D,T}(ntuple(_ -> rand(task_rng, T) - T(0.5), D))
                            cm = SVector{D,T}(uc * frac) + center
                        else
                            cm = _random_cm_in_cell(task_rng, cm_lo, cm_extent, has_pbc, uc, center, wigner_seitz)
                        end
                        angles = _random_angles(task_rng, SVector{D,T})
                        mp = MoleculePosition(cm, angles)

                        # Never accept an overlapping trial, regardless of how
                        # good its own constraint penalty is: overlapping the
                        # fixed structure is always worse than not, and isn't
                        # reflected in constraint_penalty_sum at all.
                        if packmol_system.avoid_overlap &&
                           overlaps_fixed(mp, st.reference_coordinates, task_fixed_sys, mol_positions, fixed_lo, fixed_hi, tol)
                            continue
                        end

                        # A trial that does not satisfy the constraints is fitted to
                        # them (constraints only). The fit moves it to the nearest
                        # point of its region, so taking every fitted trial would
                        # concentrate the molecules where the regions are nearest to
                        # the most random points (e.g. the directions of the corners
                        # of the box of the bounds, for a spherical shell). A fitted
                        # trial is accepted only if its center of mass moved less than
                        # `max_shift`: a trial near the region, accepted with the same
                        # probability everywhere along it, as a uniform sample of it.
                        fx = constraint_penalty_sum(mp, st, packmol_system)
                        shift = zero(T)
                        if fx >= ok_precision && fit && itry > max_guess_try
                            fit_mp, fit_fx = _fit_to_constraints!(
                                task_x, task_mol_list, packmol_system, imol, mp, task_fg_outputs[ichunk],
                                task_atom_positions[ichunk], mol_structure_type, mol_iat_first, fit_precision, opt_nit,
                                fit_bounds,
                            )
                            if fit_fx < fx && !(packmol_system.avoid_overlap &&
                                overlaps_fixed(fit_mp, st.reference_coordinates, task_fixed_sys, mol_positions, fixed_lo, fixed_hi, tol))
                                shift = norm(fit_mp.cm - mp.cm)
                                mp, fx = fit_mp, fit_fx
                            end
                        end
                        key = fx < ok_precision ? (0, shift) : (1, fx)
                        if key < best_key
                            best_key = key
                            best_mp = mp
                        end
                        if fx < ok_precision
                            if shift <= max_shift
                                task_accepted[ichunk] += 1
                                break
                            end
                        end
                    end
                    best_key[1] == 0 && (task_feasible[ichunk] += 1)
                    packmol_system.molecule_positions[imol] = best_mp
                    next!(progress_meter)
                end
            end
        end
        imol_offset += st.number_of_molecules
    end
    finish!(progress_meter)
    _println("  Trials accepted: ", sum(task_accepted), "/", nfree, " (satisfying the constraints: ", sum(task_feasible), ")")
    return packmol_system
end

#
# Fit the trial position `mp` of molecule `imol` to its constraints only (as
# `_movebad_place_molecule!`), with the buffers `x` (2D values) and `mol_list` (one
# element) of the caller. Returns the fitted position and its constraint penalty.
#
function _fit_to_constraints!(
    x::Vector{T}, mol_list::Vector{Int}, packmol_system::PackmolSystem{D,T}, imol::Int,
    mp::MoleculePosition{D,T}, fg_output::InteratomicDistanceFG{D,T}, atom_positions::Vector{SVector{D,T}},
    mol_structure_type::Vector{Int}, mol_iat_first::Vector{Int}, precision::T, opt_nit::Int,
    bounds::Tuple{Vector{T},Vector{T}},
) where {D,T}
    mol_list[1] = imol
    reinterpret(MoleculePosition{D,T}, x)[1] = mp
    lower, upper = bounds
    spgresult = spgbox!(
        (g, x) -> constraint_only_fg_for_mols!(
            g, x, fg_output, packmol_system, atom_positions, mol_list, mol_structure_type, mol_iat_first,
        ),
        x;
        lower, upper,
        nitmax=opt_nit,
        nfevalmax=10 * opt_nit,
        callback=(result) -> result.f < precision,
    )
    return reinterpret(MoleculePosition{D,T}, x)[1], T(spgresult.f)
end

#
# One overlap-check system per task: `fixed_sys` itself, and new ones built
# like it for the other tasks (rather than deep copies of it, because
# `deepcopy` does not compile with `juliac --trim`).
#
_task_overlap_check_systems(::Nothing, packmol_system, tol, ntasks) = nothing
function _task_overlap_check_systems(fixed_sys::S, packmol_system, tol, ntasks) where {S}
    task_fixed_syss = S[fixed_sys]
    for _ in 2:ntasks
        new_sys = first(_build_overlap_check_system(packmol_system, tol))
        new_sys isa S && push!(task_fixed_syss, new_sys)
    end
    return task_fixed_syss
end

#
# The bounds of the centers of mass of a structure type (`compute_cm_bounds`) are
# those of its molecules, which the constraint fits leave on the near side of the
# region of the constraints, and rarely at its extremes along the axes. For a thin
# region (as a spherical shell of rigid lipids), the faces of that box cut it, and
# the trials accepted near the region (see `_place_molecules!` and
# `_movebad_place_molecule!`) would be rarer, or absent, there. The random trials are
# drawn in the bounds extended by the radius of the molecule, within the periodic cell.
#
function _extended_cm_bounds(
    packmol_system::PackmolSystem{D,T}, st::StructureType{D,T}, lo::SVector{D,T}, hi::SVector{D,T},
) where {D,T}
    isempty(st.reference_coordinates) && return lo, hi
    r = maximum(norm, st.reference_coordinates)
    lo, hi = lo .- r, hi .+ r
    if !isnothing(packmol_system.unitcell)
        half = _cm_half_extent(packmol_system)
        center = packmol_system.unitcell_center
        lo, hi = max.(lo, center - half), min.(hi, center + half)
    end
    return lo, hi
end
