"""
    packmol(input_file::String; kargs...)

Read a Packmol input file, run the packing optimization, and write the output.
Returns `true` if the packing converged within `nloop` loops, `false` otherwise.
"""
function packmol(input_file::String; D::Int=3, T::DataType=Float64, kargs...)
    packmol_system = read_packmol_input(input_file; D, T)
    return packmol(packmol_system; kargs...)
end

"""
    packmol(packmol_system::PackmolSystem; kargs...)

Run the packing optimization on a `PackmolSystem`.
Returns `true` if the packing converged within `nloop` loops, `false` otherwise.
"""
function packmol(
    packmol_system::PackmolSystem{D,T};
    parallel::Bool=true,
    iprint::Int=10,
    nloop::Int=200,
    maxit::Union{Nothing,Int}=nothing,
    movefrac::T=T(0.05),
    n_stall_iterations::Int=10,
    f_stall_tolerance::T=T(0.10),
    seed::Int=packmol_system.seed,
    restart::Bool=false,
) where {D,T}
    maxit = something(maxit, 800)
    # Initialize RNG and molecule positions. Negative seed (Fortran Packmol's
    # `seed -1`, used by the Recipes as their default) means "pick a random
    # seed" rather than a literal RNG seed — passing it straight to
    # Xoshiro(seed) throws a DomainError on Julia 1.10 (LTS), since negative
    # integers there hit make_seed's "n must be non-negative" path.
    RNG = seed < 0 ? Random.Xoshiro() : Random.Xoshiro(seed)

    tstart = time()    

    # Print title
    _version = pkgversion(Packmol)
    println()
    println(hash_line)
    println("  PACKMOL - Packing optimization for the automated generation of")
    println("  starting configurations for molecular dynamics simulations.")
    println()
    @printf("%62s\n", "Version $_version ")
    println(hash_line)
    println()

    # Build index of free (non-fixed) molecules
    free_mol_indices = Int[]
    imol = 0
    for st in packmol_system.structure_types
        for _ in 1:st.number_of_molecules
            imol += 1
            if !st.fixed.fixed
                push!(free_mol_indices, imol)
            end
        end
    end
    nfree = length(free_mol_indices)

    # Pre-allocate scratch buffers once; reused across all hot-path calls.
    buffers = MemoryBuffers(packmol_system)

    # Restart (Fortran's `restart_from`): read saved molecule positions
    # instead of deriving a starting point the normal way, for either the
    # whole system or specific structure types — fixed molecules are never
    # included, since their position is already fully determined by their
    # own `fixed` keyword regardless of restart. A whole-system restart_from
    # skips the initial-approximation pipeline entirely for the free
    # molecules it covers (that's the point — for a large system, it's the
    # expensive part); a per-structure-type one is applied on top, after
    # that pipeline runs (whether or not it was itself skipped), since
    # re-placing only some molecule types still needs the others normally
    # initialized. Reference coordinates must already be centered (matching
    # what `_align_molecule` assumes) before either is read, since a
    # PDB-based restart_from needs them for its rigid-body alignment.
    _center_reference_coordinates!(packmol_system)
    if !isnothing(packmol_system.restart_from)
        println("  Restarting all free molecules from: ", packmol_system.restart_from)
        segments = Tuple{Int,Int,Vector{SVector{D,T}}}[
            (st.number_of_molecules, st.natoms, st.reference_coordinates)
            for st in packmol_system.structure_types if !st.fixed.fixed
        ]
        positions = _restart_positions(packmol_system.restart_from, segments, MoleculePosition{D,T})
        for (k, imol) in enumerate(free_mol_indices)
            packmol_system.molecule_positions[imol] = positions[k]
        end
        restart = true
    end

    # Pre-optimization: set initial approximation (random placement + constraint fitting)
    set_initial_approximation!(packmol_system, free_mol_indices, RNG; restart, buffers)

    imol_offset = 0
    for st in packmol_system.structure_types
        if !st.fixed.fixed && !isnothing(st.restart_from)
            println("  Restarting structure type '", basename(st.filename), "' from: ", st.restart_from)
            segments = [(st.number_of_molecules, st.natoms, st.reference_coordinates)]
            positions = _restart_positions(st.restart_from, segments, MoleculePosition{D,T})
            packmol_system.molecule_positions[imol_offset+1:imol_offset+st.number_of_molecules] .= positions
        end
        imol_offset += st.number_of_molecules
    end

    # check mode: write the initial approximation and return
    if packmol_system.check
        println("  Check mode: writing initial approximation to $(packmol_system.output_file)")
        if !isempty(packmol_system.output_file)
            write_output(packmol_system)
        end
        return false
    end

    # Compute initial atom positions (reuse pre-allocated buffer)
    natoms = length(packmol_system.atoms)
    atom_positions = buffers.atom_positions
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)

    # Determine unit cell for CellListMap
    has_pbc = !isnothing(packmol_system.unitcell)
    if has_pbc
        # PBC mode: use the actual unit cell
        unitcell = packmol_system.unitcell
    else
        # Non-PBC mode: inflate bounding box so CellListMap treats it as a large box.
        # Use constraint-derived CM bounds when available, to avoid a huge unitcell
        # when initial constraint adjustment didn't fully converge (some molecules
        # may still be at sidemax coordinates).
        cm_lo, cm_hi = compute_cm_bounds(packmol_system)
        # Overall bounds across all structure types
        all_lo = reduce((a,b) -> min.(a,b), cm_lo)
        all_hi = reduce((a,b) -> max.(a,b), cm_hi)
        # Add margin for molecule extent (max radius of reference coordinates)
        max_extent = zero(T)
        for st in packmol_system.structure_types
            for r in st.reference_coordinates
                max_extent = max(max_extent, norm(r))
            end
        end
        margin = max_extent + packmol_system.tolerance
        lo = all_lo .- margin
        hi = all_hi .+ margin
        box_size = hi - lo
        # Guard: if bounds are degenerate, fall back to bounding box of atom positions
        if any(box_size .≤ zero(T))
            lo, hi = compute_bounding_box(atom_positions)
            box_size = hi - lo
        end
        unitcell = T(1.2) * box_size
    end

    # CellListMap cutoff: at least the packing tolerance,
    # adjusted so the number of cells doesn't vastly exceed the number of atoms
    tol = packmol_system.tolerance
    packing_tol = tol + tol / 10
    if unitcell isa AbstractMatrix
        volume = abs(det(unitcell))
    else
        volume = prod(unitcell)
    end
    cutoff = _capped_cutoff(volume, packing_tol, natoms, D)

    # Set up CellListMap
    fg_output = InteratomicDistanceFG{D,T}(packmol_system)
    println("  Total number of atoms: ", natoms)
    println("  Number of free molecules: ", nfree)
    println("  Number of variables: ", nfree * 2 * D)
    cl_system = ParticleSystem(
        xpositions=atom_positions,
        unitcell=unitcell,
        cutoff=cutoff,
        output=fg_output,
        output_name=:fg,
        parallel=parallel,
    )

    # Set up optimization variables: only free molecule DOFs (reuse pre-allocated buffer)
    x = buffers.x
    x_mol = reinterpret(MoleculePosition{D,T}, x)
    for (k, imol) in enumerate(free_mol_indices)
        x_mol[k] = packmol_system.molecule_positions[imol]
    end
    auxvecs = buffers.vaux

    # Placement region for movebad! randomization
    mol_structure_type = _build_mol_structure_type(packmol_system)
    precision = packmol_system.tolerance_precision

    # `constrain_rotation` bounds (per structure type, per axis) as hard
    # bounds on the rotation-angle optimization variables — matching Fortran
    # Packmol's own `pgencan`, which sets GENCAN's `l`/`u` this way rather
    # than adding a soft penalty term. Built once (bounds don't change loop
    # to loop) in the same flat MoleculePosition-reinterpreted layout as `x`;
    # translation DOFs are always unbounded. `nothing` (not just ±Inf
    # vectors) when no structure type constrains any axis, so SPGBox skips
    # the bound-checking overhead entirely in the common case.
    any_rotation_constrained = any(packmol_system.structure_types) do st
        any(!isnothing, st.rotation_bounds)
    end
    lower, upper = if any_rotation_constrained
        lower_mol = Vector{MoleculePosition{D,T}}(undef, nfree)
        upper_mol = Vector{MoleculePosition{D,T}}(undef, nfree)
        cm_lo = SVector{D,T}(ntuple(_ -> T(-Inf), D))
        cm_hi = SVector{D,T}(ntuple(_ -> T(Inf), D))
        for (k, imol) in enumerate(free_mol_indices)
            bounds = packmol_system.structure_types[mol_structure_type[imol]].rotation_bounds
            ang_lo = SVector{D,T}(ntuple(d -> isnothing(bounds[d]) ? T(-Inf) : bounds[d][1], D))
            ang_hi = SVector{D,T}(ntuple(d -> isnothing(bounds[d]) ? T(Inf) : bounds[d][2], D))
            lower_mol[k] = MoleculePosition(cm_lo, ang_lo)
            upper_mol[k] = MoleculePosition(cm_hi, ang_hi)
        end
        reinterpret(T, lower_mol), reinterpret(T, upper_mol)
    else
        nothing, nothing
    end

    # Outer packing loop (following Fortran Packmol gencanloop):
    # Each iteration runs a short optimization, evaluates per-molecule
    # contributions, and randomly re-places the worst molecules.
    println()
    println(dash_line)
    println("  Packing $nfree free molecules ($(packmol_system.nmols) total)...")
    println(dash_line)
    # Evaluate and print initial function value
    g0 = similar(x)
    # Start packing with a looser-than-required tolerance (matches the
    # original Fortran Packmol's `discale` heuristic, default 1.1): atom
    # radii are effectively inflated by radscale in the optimization
    # target, giving the optimizer an easier target while far from
    # feasible. radscale decays toward 1.0 as improvement stalls (see the
    # end of the loop body below); true convergence checks (tol_ok, dmin)
    # are unaffected since they measure the real, unscaled distances.
    radscale = packmol_system.radscale
    f0 = fg!(g0, x, cl_system, packmol_system, atom_positions, free_mol_indices, radscale)
    @printf("  Objective function at initial point: %10.5e\n", f0)
    bestf = typemax(T)
    flast = typemax(T)
    converged = false
    best_positions = copy(packmol_system.molecule_positions)
    for loop in 0:nloop

        println()
        println(dash_line)
        @printf("  Starting packing loop: %8d\n", loop)
        @printf("  Tolerance in this loop: %8.4f\n", radscale * tol)
        println()

        # The block below runs one SPGBox chunk and classifies how it ended.
        # Most of the time that's the end of this iteration of the visible
        # `for loop` — but when a chunk ends early (packmol_callback's own
        # plateau flags, or SPGBox's internal convergence) *and* the
        # whole-chunk fimprov shows that was a false alarm (no actionable
        # stall, real progress happening), retrying immediately from the same
        # point is strictly better than falling through: printing a full
        # "Packing loop ended"/stats block and opening a new visible loop for
        # what is, from the optimization's point of view, an uninterrupted
        # continuation. So this silently loops (fresh SPGBox call, same x,
        # same radscale, no movebad!/radscale-decay, no `loop` incrementing)
        # until it hits a chunk that's either genuinely done (converged),
        # genuinely stalled (an action is warranted), or ran its full
        # maxit/nfevalmax budget — only then falling through to the reporting
        # and decision logic below. Capped at `max_silent_retries` purely as
        # a safety net against a pathological run of ever-smaller-but-still->
        # threshold gains that never falls through on its own.
        max_silent_retries = 20
        local optresult, fx, dmin, max_const, tol_ok, const_ok, fimprov, fimp_last, improved
        local stall_reasons, ambiguous_stall, raw_dmin_stall, raw_constraint_stall
        local dmin_stalled, constraint_stalled, chunk_stalled
        silent_retries = 0
        while true
            # Run a short optimization (maxit iterations per loop)
            progress_meter = Progress(maxit; desc=" Iterations: ", barlen=47)
            fg_closure = (g, x) -> fg!(g, x, cl_system, packmol_system, atom_positions, free_mol_indices, radscale)
            # Fresh stall detectors per chunk: plateau state from the *previous*
            # chunk isn't meaningful here, since movebad! (or the initial
            # placement) just changed the starting point. `dmin_stalled_flag`/
            # `constraint_stalled_flag` are set by `packmol_callback` when it
            # cuts the chunk short because of a stall in that metric (as opposed
            # to full convergence or simply exhausting the maxit budget) — the
            # outer loop below uses them, instead of an inter-loop function-value
            # comparison, both to decide whether to move bad molecules or
            # loosen/tighten radscale, and to report why this chunk ended.
            dmin_stall_detector = StallDetector{T}(n_stall_iterations)
            constraint_stall_detector = StallDetector{T}(n_stall_iterations)
            # dmin/max_const are each a worst-case-over-all-atoms extremum, not
            # an aggregate, so whichever single pair currently holds that worst
            # value routinely sits frozen for a while — stuck behind other
            # molecules resolving their own, larger violations first — even
            # while the coupled objective f is still falling sharply overall.
            # That made every chunk on a large, realistic system report "stalled"
            # within ~n_stall_iterations regardless of how much real progress f
            # was making elsewhere. `f_stall_detector` tracks f over a much
            # longer window (maxit/20, e.g. 40 iterations at the default
            # maxit=800) — long enough to average out that per-pair noise — and
            # packmol_callback only honors a dmin/const stall once f *also* shows
            # no net improvement (< f_stall_tolerance, e.g. 10%) over that
            # window: genuine overall progress vetoes a worst-case-metric stall.
            f_stall_detector = StallDetector{T}(max(1, maxit ÷ 20))
            dmin_stalled_flag = Ref(false)
            constraint_stalled_flag = Ref(false)
            # f as of the first callback invocation in this chunk (typemax(T)
            # sentinel until then) — the reference packmol_callback compares
            # against to see whether the *whole chunk so far*, not just its
            # trailing f_stall_detector window, has made significant progress.
            # See the comment above f_not_improving in packmol_callback for why
            # the narrower window alone isn't enough to justify cutting a chunk
            # short.
            f_chunk_start = Ref(typemax(T))
            optresult = spgbox!(
                fg_closure,
                x;
                callback=(result) -> packmol_callback(cl_system, tol, iprint,
                    packmol_system.tolerance_precision, packmol_system.constraint_precision, progress_meter;
                    dmin_stall_detector, constraint_stall_detector, stall_tolerance=packmol_system.stall_tolerance,
                    f_stall_detector, f=result.f, f_stall_tolerance,
                    dmin_stalled_flag, constraint_stalled_flag, f_chunk_start,
                ),
                vaux=auxvecs,
                nitmax=maxit,
                nfevalmax=10 * maxit,
                lower, upper,
            )

            # Update molecule positions with optimized values
            x_mol = reinterpret(MoleculePosition{D,T}, x)
            for (k, imol) in enumerate(free_mol_indices)
                packmol_system.molecule_positions[imol] = x_mol[k]
            end

            # Statistics: compute improvement of this loop relative to flast
            fx = optresult.f
            dmin = min(cl_system.fg.dmin, cl_system.cutoff)
            fimprov = bestf < typemax(T) ? clamp(-100 * (fx - bestf) / bestf, T(-99.99), T(99.99)) : T(100)
            # On the very first loop, flast is still its typemax sentinel (no
            # previous loop to compare against): fall back to the improvement
            # from the best function value instead of computing a bogus
            # Inf/Inf = NaN ratio against the sentinel.
            fimp_last = flast < typemax(T) ? clamp(-100 * (fx - flast) / flast, T(-99.99), T(99.99)) : fimprov
            improved = fx < bestf
            if improved
                bestf = fx
                copyto!(best_positions, packmol_system.molecule_positions)
            end
            flast = fx
            max_const = cl_system.fg.max_constraint_penalty

            # Check convergence: both tolerance and constraint precisions must be satisfied
            tol_ok = tol - dmin < packmol_system.tolerance_precision
            const_ok = max_const < packmol_system.constraint_precision

            # Report why this chunk's optimization ended: fully converged, one
            # (or more) of the stall conditions cut it short, or the chunk simply
            # ran through its maxit/nfevalmax budget. A stall condition fires
            # from either of the two dmin/constraint plateau detectors above, or
            # from SPGBox reaching its own internal convergence criterion (small
            # projected gradient) short of Packmol's own tolerance/constraint
            # targets — this happens when re-optimizing from the exact same
            # starting point produces so few iterations that the plateau
            # detectors never get the several consecutive readings they need to
            # fire, even though nothing is actually still improving. In that
            # ambiguous case, attribute it to whichever of tol_ok/const_ok is
            # still unmet (possibly both), so it feeds movebad!/radscale decay
            # below the same way an explicit plateau detection would.
            stall_reasons = String[]
            dmin_stalled_flag[] && push!(stall_reasons, "minimum distance plateaued")
            constraint_stalled_flag[] && push!(stall_reasons, "constraint violation plateaued")
            optimizer_converged_internally = optresult.ierr == 0
            # `ierr == 0` alone doesn't distinguish genuine internal convergence
            # from a callback-triggered early stop: SPGBox returns the very same
            # ierr==0 when the callback returns true (spgbox_main.jl's `return
            # SPGBoxResult(...,0,true)` on a stall-detector hit) as when it
            # actually converges. So `ambiguous_stall` — the fallback below —
            # must only apply when *neither* plateau detector already explains
            # why the chunk stopped; otherwise a pure dmin-only stall (with
            # constraints still unmet, as they almost always are while dmin is
            # the one being fixed) would get relabeled as also constraint-stalled
            # just because ierr==0 and const_ok happens to still be false.
            ambiguous_stall = optimizer_converged_internally && isempty(stall_reasons)
            if ambiguous_stall
                !tol_ok && push!(stall_reasons, "optimizer converged internally short of the distance tolerance")
                !const_ok && push!(stall_reasons, "optimizer converged internally short of the constraints")
            end
            # dmin_stalled_flag/constraint_stalled_flag (and the ambiguous_stall
            # fallback) are set from packmol_callback's *internal* trailing-window
            # view of this chunk (its own f_stall_detector, sized to maxit/20 —
            # e.g. the last ~40 iterations of an 800-iteration chunk). That's the
            # right granularity for "stop burning iterations on a plateaued
            # metric", but it's blind to the chunk as a whole: a chunk that fell
            # sharply in its first two-thirds and only flattened out in that
            # trailing window has made excellent overall progress, even though
            # the window itself looks stalled. Concretely: dmin/max_const
            # plateaued while fimprov (the real, whole-chunk improvement vs. the
            # best f seen so far, computed above) was ~99.78% — nothing here
            # actually needs fixing, the chunk should just be given another
            # normal iteration budget. So a raw plateau only becomes an
            # actionable stall (worth radscale decay or relocating molecules)
            # once the whole-chunk fimprov *also* fails to clear
            # f_stall_tolerance — genuine chunk-level progress vetoes a
            # trailing-window plateau, mirroring how f_stall_detector already
            # vetoes dmin/const plateaus with its own (shorter) window.
            #
            # This must compare against bestf (fimprov), not against the
            # immediately preceding loop's f (fimp_last): after a movebad! has
            # scrambled positions, the next chunk's f is often still far worse
            # than bestf even once it has clawed back >f_stall_tolerance of that
            # self-inflicted damage. fimp_last would then read as "significant
            # progress" purely from recovering off a low bar it created itself,
            # vetoing further movebad!/radscale action while the packing remains
            # deeply regressed relative to its actual best-ever state — visible
            # as a "stalled" -> "plateau ignored" -> "stalled" thrash with no net
            # improvement. fimprov, anchored to bestf, isn't corrupted by that
            # recent self-inflicted dip.
            significant_progress = fimprov > T(100) * f_stall_tolerance
            raw_dmin_stall = dmin_stalled_flag[] || (ambiguous_stall && !tol_ok)
            raw_constraint_stall = constraint_stalled_flag[] || (ambiguous_stall && !const_ok)
            dmin_stalled = raw_dmin_stall && !significant_progress
            constraint_stalled = raw_constraint_stall && !significant_progress
            chunk_stalled = dmin_stalled || constraint_stalled
            finish!(progress_meter)

            # No actionable stall, but *something* did cut this chunk short
            # (raw_dmin_stall/raw_constraint_stall): retry silently rather than
            # reporting and opening a new visible loop, per the reasoning above.
            no_action_plateau = !(tol_ok && const_ok) && !chunk_stalled &&
                (raw_dmin_stall || raw_constraint_stall)
            if no_action_plateau && silent_retries < max_silent_retries
                silent_retries += 1
                continue
            end
            break
        end
        loop_end_reason = if tol_ok && const_ok
            "converged"
        elseif chunk_stalled
            "stalled (" * join(stall_reasons, ", ") * ")"
        elseif raw_dmin_stall || raw_constraint_stall
            "plateau ignored (f still improving)"
        elseif optresult.ierr == 2
            "chunk function-evaluation budget (nfevalmax) reached"
        else
            "chunk iteration budget (maxit) reached"
        end

        @printf("\n  Packing loop ended: %s\n", loop_end_reason)
        @printf("  Function value from last loop: f = %10.5e\n", fx)
        @printf("  Best function value before: f = %10.5e\n", bestf)
        @printf("  Improvement from best function value: %8.2f %%\n", fimprov)
        @printf("  Improvement from last loop: %8.2f %%\n", fimp_last)
        @printf("  Maximum violation of target distance: %12.6f\n", tol - dmin)
        @printf("  Maximum violation of the constraints: %10.5e\n", max_const)

        if tol_ok && const_ok
            println()
            println(hash_line)
            @printf("\n%s Success! \n", " "^32)
            @printf("%s Final objective function value: %10.5e\n", " "^13, fx)
            @printf("%s Maximum violation of target distance: %10.6f\n", " "^13, tol - dmin)
            @printf("%s Maximum violation of the constraints: %10.5e\n", " "^13, max_const)
            println()
            println(dash_line)
            println()
            println("$(repeat(" ", 13)) Please cite this work if Packmol was useful: ")
            println()
            println("$(repeat(" ", 10))  L. Martinez, R. Andrade, E. G. Birgin, J. M. Martinez, ")
            println("$(repeat(" ", 8))  PACKMOL: A package for building initial configurations for")
            println("$(repeat(" ", 18)) molecular dynamics simulations. ")
            println("$(repeat(" ", 7))  Journal of Computational Chemistry, 30(13) pp. 2157-2164, 2009.")
            println("$(repeat(" ", 17)) https://doi.org/10.1002/jcc.21224")
            println()
            println(hash_line)
            converged = true
            break
        end

        println()
        println(dash_line)

        # Write best solution so far to output file
        if improved && !isempty(packmol_system.output_file)
            saved_positions = packmol_system.molecule_positions
            packmol_system.molecule_positions = best_positions
            write_output(packmol_system)
            packmol_system.molecule_positions = saved_positions
            println()
            println("  Current solution written to file: ", packmol_system.output_file)
        end

        # Two stall responses, mutually exclusive per loop:
        #  - Distance-related stall while radscale still has room to relax
        #    (> 1.0): ease the pairwise-distance target further rather than
        #    relocating molecules — a molecule that's merely a bit too close
        #    under the current (inflated) tolerance may resolve on its own
        #    once that tolerance is loosened, so this is tried first.
        #  - Everything else that's still stalled once radscale is already
        #    fully tightened (== 1.0), or a constraint-region stall (no
        #    amount of radscale/tolerance easing frees a molecule stuck
        #    outside its assigned region — relocating it is the only fix):
        #    move the offending molecules to new random positions.
        # `dmin_stalled`/`constraint_stalled` already require the whole-chunk
        # function improvement to also be small (see above), so acting on
        # them here can't disrupt a loop that's genuinely still improving.
        do_radscale_decay = dmin_stalled && radscale > one(T)
        do_movebad = chunk_stalled && !do_radscale_decay
        if do_radscale_decay
            radscale = max(T(0.9) * radscale, one(T))
        elseif do_movebad
            cm_min, cm_max = compute_cm_bounds(packmol_system)
            nmoved = movebad!(
                packmol_system, cl_system.fg.fmol, free_mol_indices, mol_structure_type, RNG;
                movefrac, precision,
                cm_lo_type=cm_min, cm_hi_type=cm_max,
            )
            if nmoved > 0
                println("  Moved $nmoved bad molecules randomly to new positions.")
            end
        end
        # Re-pack optimizer variables from (possibly moved) molecule positions
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(free_mol_indices)
            x_mol[k] = packmol_system.molecule_positions[imol]
        end
    end

    if !converged
        println()
        println(hash_line)
        @printf("  WARNING: packing did not converge after %d loops (best f = %.4e)\n", nloop, bestf)
        println(hash_line)
    end

    # Restore best molecule positions
    copyto!(packmol_system.molecule_positions, best_positions)

    # For PBC: wrap atom positions to the unit cell centered at unitcell_center
    if has_pbc
        center = packmol_system.unitcell_center
        compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)
        for i in eachindex(atom_positions)
            atom_positions[i] = wrap_to_center(atom_positions[i], packmol_system.unitcell, center)
        end
        # Update molecule centers of mass to wrapped positions
        # (recompute CM from the wrapped atom positions for each molecule)
        iat = 0
        imol = 0
        for st in packmol_system.structure_types
            for _ in 1:st.number_of_molecules
                imol += 1
                mol_cm = zero(SVector{D,T})
                for j in 1:st.natoms
                    iat += 1
                    mol_cm += atom_positions[iat]
                end
                mol_cm /= st.natoms
                packmol_system.molecule_positions[imol] = MoleculePosition(
                    mol_cm, packmol_system.molecule_positions[imol].angles
                )
            end
        end
    end

    # Write output file if specified
    if !isempty(packmol_system.output_file)
        write_output(packmol_system)
        println()
        println("  Solution written to file: ", packmol_system.output_file)
    end

    println()
    println(dash_line)
    tend = time()
    @printf("  Running time: %12.4f seconds.\n", tend - tstart)
    println(dash_line)
    println()

    return converged
end

#
# Aesthetic line constants (matching Fortran Packmol)
#
const dash_line = repeat('-', 80)
const hash_line = repeat('#', 80)

#
# Tracks whether a scalar convergence metric has stopped meaningfully
# improving over a trailing *window* of SPGBox iterations within a single
# packing-loop chunk. Used for two independent metrics: the minimum
# interatomic distance `dmin` (the quantity `tol_ok` depends on, larger is
# better) and the maximum constraint violation `max_constraint_penalty` (the
# quantity `const_ok` depends on, smaller is better). A handful of molecules
# genuinely stuck (deep inside a fixed structure, or wedged against others
# with nowhere left to go) show up as one of these flatlining, even while
# the coupled many-body objective f may still be improving elsewhere in the
# system (tried f itself and reverted: on this coupled objective, an
# instantaneous gradient/function ratio can look artificially small for the
# first iteration or two after a chunk starts, e.g. right after movebad!
# scrambles some positions, well before the spectral step-size estimate has
# calibrated — that caused chunks to bail out almost immediately, before
# making any real progress).
#
# A fixed-size circular buffer holds the last `n_stall_iterations` values.
# `is_stalled!` compares the current value against the one from exactly
# `n_stall_iterations` calls ago and reports true only once that *net* change
# over the whole window falls below `rel_tol` — not once any single step
# does. A metric that improves slowly but steadily (each individual SPGBox
# step under `rel_tol`, e.g. because the tolerance is tight relative to the
# step size) still shows real progress once accumulated over the window, and
# is correctly not flagged as stalled; a metric genuinely flatlined shows
# ~zero net change over the same window regardless of how it's chopped up.
# Before the buffer has seen `n_stall_iterations` values there's no full
# window to compare yet, so the metric is never counted as stalled.
#
mutable struct StallDetector{T}
    window::Vector{T}
    head::Int  # next slot to write (1-based); also the slot holding the oldest value once full
    count::Int  # total values seen so far, capped at length(window)
    StallDetector{T}(n_stall_iterations::Int) where {T} = new{T}(Vector{T}(undef, n_stall_iterations), 1, 0)
end

function is_stalled!(
    detector::StallDetector{T}, value::T;
    rel_tol::T, larger_is_better::Bool,
) where {T}
    n = length(detector.window)
    detector.count += 1
    if detector.count <= n
        detector.window[detector.head] = value
        detector.head = detector.head == n ? 1 : detector.head + 1
        return false
    end
    oldest = detector.window[detector.head]
    detector.window[detector.head] = value
    detector.head = detector.head == n ? 1 : detector.head + 1
    rel_improvement = if oldest > zero(T)
        larger_is_better ? (value - oldest) / oldest : (oldest - value) / oldest
    else
        zero(T)
    end
    return rel_improvement < rel_tol
end

#
# SPGBox callback: print progress and check convergence
#
function packmol_callback(
    cl_system, tol, iprint, tolerance_precision, constraint_precision, progress_meter;
    dmin_stall_detector::Union{Nothing,StallDetector} = nothing,
    constraint_stall_detector::Union{Nothing,StallDetector} = nothing,
    stall_tolerance = 1e-2 * tolerance_precision,
    f_stall_detector::Union{Nothing,StallDetector} = nothing,
    f::Union{Nothing,Real} = nothing,
    f_stall_tolerance = 0.10,
    dmin_stalled_flag::Union{Nothing,Ref{Bool}} = nothing,
    constraint_stalled_flag::Union{Nothing,Ref{Bool}} = nothing,
    f_chunk_start::Union{Nothing,Ref} = nothing,
)
    dmin = min(cl_system.fg.dmin, cl_system.cutoff)
    max_const = cl_system.fg.max_constraint_penalty
    next!(progress_meter; showvalues = [
        (" Minimum distance: ", dmin),
        (" Maximum constraint violation: ", max_const),
    ])
    tol_ok = tol - dmin < tolerance_precision
    const_ok = max_const < constraint_precision
    if tol_ok && const_ok
        return true
    end
    # dmin/max_const are each a worst-case-over-all-atoms extremum, not an
    # aggregate: in a large system it's common for whichever single pair
    # currently holds that worst value to sit frozen for a stretch — stuck
    # behind other molecules resolving their own, larger violations first —
    # even while the coupled objective f is still falling sharply overall.
    # So neither metric's own plateau is trusted on its own: it only counts
    # as a real stall once f *also* shows no significant net improvement
    # (< f_stall_tolerance, e.g. 10%) over its own, much longer window
    # (`f_stall_detector`, sized to maxit/20 by the caller) — genuine
    # overall progress vetoes a worst-case-metric plateau. `is_stalled!` is
    # still called unconditionally each iteration (not short-circuited) so
    # every detector's window stays populated regardless of which branch
    # ends up mattering.
    f_window_stalled = isnothing(f_stall_detector) || isnothing(f) ||
        is_stalled!(f_stall_detector, oftype(dmin, f); rel_tol=oftype(dmin, f_stall_tolerance), larger_is_better=false)
    # f_stall_detector's window (~40 iterations at the default maxit=800) is
    # sized to average out per-pair noise, but it's still only a trailing
    # slice of the chunk: a chunk that fell sharply over its first several
    # hundred iterations and has merely leveled off in this recent window has
    # made excellent progress overall, and cutting it short here would just
    # force the outer loop to close this chunk out and immediately open a new
    # one from the same point — wasted bookkeeping for no benefit, since nitmax
    # iterations remain unused. So a plateaued trailing window is itself
    # vetoed by whole-chunk progress: `f_chunk_start` is set once, from the
    # first f this callback ever sees in this chunk, and as long as f has
    # dropped more than `f_stall_tolerance` from that opening value, the chunk
    # keeps running regardless of what the trailing window shows.
    if !isnothing(f_chunk_start) && !isnothing(f) && f_chunk_start[] == typemax(f_chunk_start[])
        f_chunk_start[] = oftype(dmin, f)
    end
    chunk_progress_significant = !isnothing(f_chunk_start) && !isnothing(f) &&
        f_chunk_start[] < typemax(f_chunk_start[]) && f_chunk_start[] > zero(f_chunk_start[]) &&
        (f_chunk_start[] - oftype(dmin, f)) / f_chunk_start[] > oftype(dmin, f_stall_tolerance)
    f_not_improving = f_window_stalled && !chunk_progress_significant
    # Cut this chunk short once either the minimum distance or the maximum
    # constraint violation shows ~no net improvement over the trailing
    # `n_stall_iterations`-iteration window (and, per above, f isn't picking
    # up the slack) — grinding through the rest of this chunk's maxit budget
    # on a stuck molecule cannot help; movebad! (once this chunk returns to
    # the outer loop) is what actually relocates them. Each metric's stall
    # detector keeps tracking regardless, but only gates the cut-short
    # decision (and the corresponding flag, which the outer loop reports as
    # the reason this chunk ended) while its own criterion hasn't yet
    # converged: once dmin (or the constraints) is already within tolerance,
    # there's nothing left to stall on for it.
    dmin_stalled = !isnothing(dmin_stall_detector) &&
        is_stalled!(dmin_stall_detector, dmin; rel_tol=stall_tolerance, larger_is_better=true) &&
        f_not_improving
    const_stalled = !isnothing(constraint_stall_detector) &&
        is_stalled!(constraint_stall_detector, max_const; rel_tol=stall_tolerance, larger_is_better=false) &&
        f_not_improving
    if !tol_ok && dmin_stalled
        isnothing(dmin_stalled_flag) || (dmin_stalled_flag[] = true)
    end
    if !const_ok && const_stalled
        isnothing(constraint_stalled_flag) || (constraint_stalled_flag[] = true)
    end
    if (!tol_ok && dmin_stalled) || (!const_ok && const_stalled)
        return true
    end
    return false
end
