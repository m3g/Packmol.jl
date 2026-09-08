#=
    _radscale_at(elapsed::Int, radscale_start::T, radscale_decay_loops::T) where {T}

Loop-indexed radscale schedule: `radscale_start` at `elapsed == 0`, decaying
exponentially to `1.0` (no inflation), reaching it within `1e-3` relative
precision by `elapsed == radscale_decay_loops`, and staying at `1.0` beyond
that. `elapsed` is simply the current packing loop count (see the outer
packing loop in `packmol`) — every atom follows this same schedule, with no
per-atom or per-molecule restart.
=#
function _radscale_at(elapsed::Int, radscale_start::T, radscale_decay_loops::T) where {T}
    (radscale_start <= one(T) || radscale_decay_loops <= zero(T) || elapsed >= radscale_decay_loops) && return one(T)
    k = -log(T(1e-3)) / radscale_decay_loops
    return one(T) + (radscale_start - one(T)) * exp(-k * elapsed)
end

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
    n_stall_iterations::Int=40,
    f_stall_tolerance::T=T(0.01),
    movebad_tolerance::T=T(0.10),
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
    # whole system or specific structure types — fixed molecules' own
    # positions are never overwritten, since they're already fully
    # determined by their own `fixed` keyword regardless of restart. A
    # whole-system restart_from skips the initial-approximation pipeline
    # entirely for the free molecules it covers (that's the point — for a
    # large system, it's the expensive part); a per-structure-type one is
    # applied on top, after that pipeline runs (whether or not it was itself
    # skipped), since re-placing only some molecule types still needs the
    # others normally initialized. Reference coordinates must already be
    # centered (matching what `_align_molecule` assumes) before either is
    # read, since a PDB-based restart_from needs them for its rigid-body
    # alignment.
    #
    # A *whole-system* restart source (PDB or raw) is whatever a prior
    # `write_output`/`_write_restart_files` call wrote for the entire
    # system, so it lists every molecule of every structure type, fixed
    # ones included, interleaved in structure-declaration order — not just
    # the free ones. `segments` below therefore spans every structure type
    # (not just the non-fixed ones), with `extract` (the tuple's 4th field)
    # marking which ones to actually read a position back from; a fixed
    # segment's atoms/lines still have to be counted and skipped over so the
    # expected total and the file's actual layout line up (see
    # `_restart_positions_from_atoms`/`_restart_positions`), even though its
    # own position is discarded.
    _center_reference_coordinates!(packmol_system)
    if !isnothing(packmol_system.restart_from)
        println("  Restarting all free molecules from: ", packmol_system.restart_from)
        segments = Tuple{Int,Int,Vector{SVector{D,T}},Bool}[
            (st.number_of_molecules, st.natoms, st.reference_coordinates, !st.fixed.fixed)
            for st in packmol_system.structure_types
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
            # A per-structure-type restart source only ever contains this
            # one type's own molecules (e.g. its own `restart_to` file from
            # an earlier run), so its segment list is just itself — always
            # extracted, nothing to skip.
            segments = [(st.number_of_molecules, st.natoms, st.reference_coordinates, true)]
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
    mol_iat_first = _build_mol_iat_first(packmol_system)
    precision = packmol_system.tolerance_precision
    # Worst fmol ever observed per structure type, across every movebad! call
    # over the whole run (see movebad.jl for why this must persist rather
    # than being recomputed fresh each call).
    fmol_max_type = zeros(T, length(packmol_system.structure_types))

    # Fixed-structure overlap-check system and constraint-only scratch for
    # movebad! (see movebad.jl): built once here, since the fixed atoms never
    # move over the course of the run, and reused on every relocation.
    movebad_overlap_tol = packmol_system.radscale * packmol_system.tolerance
    movebad_fixed_sys, movebad_fixed_lo, movebad_fixed_hi =
        _build_overlap_check_system(packmol_system, movebad_overlap_tol)
    movebad_fg_output = InteratomicDistanceFG{D,T}(packmol_system)

    # `constrain_rotation` bounds (per structure type, per axis) as hard
    # bounds on the rotation-angle optimization variables — matching Fortran
    # Packmol's own `pgencan`, which sets GENCAN's `l`/`u` this way rather
    # than adding a soft penalty term. Built once (bounds don't change loop
    # to loop) in the same flat MoleculePosition-reinterpreted layout as `x`.
    #
    # Under PBC, translation is *also* hard-bounded — to exactly one
    # canonical period of the cell, centered at `unitcell_center` — rather
    # than left unbounded. Without this, a molecule can be pushed by a large
    # gradient (e.g. deep inside a non-periodic constraint's violation
    # region, like `below plane`, whose restoring force has no periodicity
    # of its own even though the cell does) past a periodic face in a single
    # SPGBox trial step. Since the wrap used to *evaluate* fg! at that trial
    # point (`wrap_to_center`) is a sawtooth function of the raw coordinate —
    # continuous almost everywhere but jumping by a full period exactly at
    # each face — a step that crosses one or more periods lands the wrapped
    # position somewhere uncorrelated with where the gradient was actually
    # aiming. That corrupts the spectral step-size estimate (built from the
    # secant/curvature ratio between consecutive gradients), which then
    # keeps producing further oversized, uncorrelated steps: confirmed by
    # instrumentation on a single stuck molecule — its *raw* cm spiralled
    # from z=49 to z=-14295 (nearly 300 box-periods away) over just a dozen
    # SPGBox iterations, each burning dozens of backtracking function
    # evaluations for near-zero net progress, before the search happened to
    # land somewhere feasible by chance. Bounding cm to one period is not a
    # loss of generality — every reachable physical configuration already has
    # a representative inside that one period, via wrapping — it just forces
    # SPGBox to find it directly instead of possibly overshooting through
    # several periods and having to claw back. `nothing` (not just ±Inf
    # vectors) when neither PBC nor any structure type's rotation is
    # constrained, so SPGBox skips the bound-checking overhead entirely in
    # that common case.
    any_rotation_constrained = any(packmol_system.structure_types) do st
        any(!isnothing, st.rotation_bounds)
    end
    lower, upper = if has_pbc || any_rotation_constrained
        lower_mol = Vector{MoleculePosition{D,T}}(undef, nfree)
        upper_mol = Vector{MoleculePosition{D,T}}(undef, nfree)
        cm_lo, cm_hi = if has_pbc
            # Axis-aligned half-extent of the (possibly triclinic) periodic
            # cell — exact for an orthorhombic (diagonal) unitcell, and a
            # (safe, if not perfectly tight) enclosing box for a sheared one.
            half_extent = SVector{D,T}(
                ntuple(i -> sum(abs(packmol_system.unitcell[i, j]) for j in 1:D) / 2, D)
            )
            packmol_system.unitcell_center .- half_extent, packmol_system.unitcell_center .+ half_extent
        else
            SVector{D,T}(ntuple(_ -> T(-Inf), D)), SVector{D,T}(ntuple(_ -> T(Inf), D))
        end
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
    # original Fortran Packmol's `discale` heuristic, default 1.1): every
    # atom's own radius is inflated by radscale in the optimization target,
    # giving the optimizer an easier target while far from feasible.
    # `atom_radii` is the mutable, per-atom working radius used by the
    # optimizer — every atom follows the same loop-indexed `_radscale_at`
    # schedule toward its own floor (`atom_radii_floor`, its user-specified
    # or default radius), based on the packing loop counter alone: a
    # molecule relocated by movebad! is not given its own separate,
    # restarted schedule — it shares whatever radscale every other atom is
    # currently at. True convergence checks (tol_ok, dmin) are unaffected
    # since they measure the real, unscaled distances, and `atom_radii_floor`
    # is what the "true objective" evaluations below (`f_true_loop_start`/
    # `f_true_loop_end`) use.
    atom_radii_floor = Vector{T}(undef, natoms)
    for (iat, a) in enumerate(packmol_system.atoms)
        atom_radii_floor[iat] = a.radius
    end
    # radscale decays from packmol_system.radscale (at loop == 0) to 1.0 by
    # loop == radscale_decay_loops == nloop, the same for every atom.
    radscale_decay_loops = T(nloop)
    atom_radii = packmol_system.radscale .* atom_radii_floor
    # Objective function at the initial point, on the same unscaled
    # (radscale == 1.0) basis as every other cross-loop comparison below.
    f0_true = fg!(g0, x, cl_system, packmol_system, atom_positions, free_mol_indices, atom_radii_floor)
    @printf("  Objective function at initial point: %10.5e\n", f0_true)
    # `bestf` starts at the true initial-point value (not a typemax
    # sentinel): loop 0's own `f_true_loop_start` is evaluated at this exact
    # same position, so `fimprov` (vs. `bestf`) and `fimp_within_loop` (vs.
    # `f_true_loop_start`) are identical at loop 0 — as they must be, since
    # no molecules have been relocated yet to make the two baselines diverge.
    # A typemax sentinel would instead force `improved = true` unconditionally
    # on loop 0 regardless of whether the true objective actually got worse,
    # silently discarding this good initial configuration in `best_positions`.
    bestf = f0_true
    converged = false
    best_positions = copy(packmol_system.molecule_positions)
    for loop in 0:nloop

        println()
        println(dash_line)
        @printf("  Starting packing loop: %8d\n", loop)
        atom_radii_lo, atom_radii_hi = extrema(atom_radii)
        @printf("  Atom radii in this loop: [ %8.4f - %8.4f ]\n", atom_radii_lo, atom_radii_hi)
        println()

        # True (floor) objective value at the position this loop starts
        # from — always evaluated at each atom's own floor radius
        # (atom_radii_floor), regardless of this loop's own (possibly
        # inflated) working atom_radii, so it can be compared on equal
        # footing with f_true_loop_end below to isolate this loop's own net
        # effect on the real objective (see "Improvement within this loop"
        # reported at the end of the loop). When atom_radii already equals
        # atom_radii_floor for every atom in this loop (fully decayed, no
        # atom freshly fattened by movebad!), the optimizer directly
        # minimizes this same quantity, so that comparison is guaranteed
        # non-increasing; otherwise, the optimizer is chasing an inflated
        # target instead, so the real objective isn't guaranteed to move in
        # either direction.
        f_true_loop_start = fg!(g0, x, cl_system, packmol_system, atom_positions, free_mol_indices, atom_radii_floor)

        # Run a short optimization (maxit iterations per loop)
        progress_meter = Progress(maxit; desc=" Iterations: ", barlen=47)
        fg_closure = (g, x) -> fg!(g, x, cl_system, packmol_system, atom_positions, free_mol_indices, atom_radii)
        # Fresh stall detectors per chunk: plateau state from the *previous*
        # chunk isn't meaningful here, since movebad! (or the initial
        # placement) just changed the starting point. `dmin_stalled_flag`/
        # `constraint_stalled_flag` are set by `packmol_callback` when it
        # cuts the chunk short because of a stall in that metric (as opposed
        # to full convergence or simply exhausting the maxit budget) — this
        # only stops this chunk's SPGBox run early so it doesn't grind
        # through the rest of its maxit/nfevalmax budget for nothing; the
        # outer loop below reports why the chunk ended from them, but decides
        # whether to move bad molecules separately, from the true
        # (floor-radii) objective's own improvement over the whole loop.
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
        # no net improvement (< f_stall_tolerance, e.g. 5%) over that
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
        # `progress_vs_true_start` (see packmol_callback) is a one-time grace
        # period, not a standing veto: this flag latches once that signal has
        # been used to excuse a plateau, so it can't keep excusing one for
        # the rest of the chunk — see the comment above it in
        # packmol_callback for why an unlatched version is a bug (a
        # one-way-true "still improving" veto that never expires).
        f_true_start_progress_used = Ref(false)
        optresult = spgbox!(
            fg_closure,
            x;
            callback=(result) -> packmol_callback(cl_system, tol, iprint,
                packmol_system.tolerance_precision, packmol_system.constraint_precision, progress_meter;
                dmin_stall_detector, constraint_stall_detector, stall_tolerance=packmol_system.stall_tolerance,
                f_stall_detector, f=result.f, f_stall_tolerance, f_true_loop_start,
                dmin_stalled_flag, constraint_stalled_flag, f_chunk_start, f_true_start_progress_used,
                nfeval=result.nfeval, nfevalmax=10 * maxit, gnorm=result.gnorm,
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

        # True (floor) objective value at the position this loop ends at,
        # paired with f_true_loop_start above: this isolates what this
        # loop's own optimization did to the real objective. Since this
        # loop's own optimizer worked under atom_radii (possibly still
        # inflated above atom_radii_floor at this loop), a negative
        # fimp_within_loop below is possible — it means the inflated target
        # it was actually chasing pulled the real, unscaled objective the
        # wrong way.
        #
        # dmin/max_const/fimprov/bestf, and the "Function value from last
        # loop" printed below, are all deliberately based on this true,
        # floor-radii value rather than the optimizer's own (possibly
        # radscale-inflated) `optresult.f`: the latter isn't on the same
        # scale from one loop to the next, since atom_radii decays toward
        # atom_radii_floor as loops progress, so a raw optresult.f-to-f0_true
        # (or optresult.f-to-bestf) comparison would compare two different
        # objective functions rather than two snapshots of the same one — and
        # printing it right next to the floor-radii "Best function value
        # before"/percentages below would be misleading for the same reason.
        # (dmin and max_const are unaffected either way — the
        # interatomic-distance minimum and the geometric constraint
        # penalties are both independent of atom_radii — so reading them off
        # this floor-radii evaluation gives exactly the same value they'd
        # have under the working atom_radii.)
        f_true_loop_end = fg!(g0, x, cl_system, packmol_system, atom_positions, free_mol_indices, atom_radii_floor)
        dmin = min(cl_system.fg.dmin, cl_system.cutoff)
        fimp_within_loop = clamp(-100 * (f_true_loop_end - f_true_loop_start) / f_true_loop_start, T(-99.99), T(99.99))
        fimprov = clamp(-100 * (f_true_loop_end - bestf) / bestf, T(-99.99), T(99.99))
        improved = f_true_loop_end < bestf
        bestf_before_loop = bestf
        if improved
            bestf = f_true_loop_end
            copyto!(best_positions, packmol_system.molecule_positions)
        end
        max_const = cl_system.fg.max_constraint_penalty
        # The call above overwrote cl_system.fg (fmol, gradients, dmin, ...)
        # with atom_radii_floor values; restore it to this loop's own working
        # atom_radii before movebad! (below) reads fmol to pick which
        # molecules to relocate.
        fg!(g0, x, cl_system, packmol_system, atom_positions, free_mol_indices, atom_radii)

        # Check convergence: both tolerance and constraint precisions must be satisfied
        tol_ok = tol - dmin < packmol_system.tolerance_precision
        const_ok = max_const < packmol_system.constraint_precision

        # Report why this chunk's optimization ended. This is purely
        # diagnostic — whether to move bad molecules is decided independently
        # below, from the true, floor-radii objective's own improvement this
        # loop (`do_movebad`), not from any of this. `dmin_stalled_flag`/
        # `constraint_stalled_flag` are set by `packmol_callback` in the
        # exact same branch where it returns `true` on a windowed dmin/const
        # plateau (having already applied its own progress veto internally,
        # via `chunk_progress_significant`) — so either flag being set is a
        # direct, reliable record that the callback itself is what ended this
        # chunk, not a post-hoc reconstruction from `optresult.ierr` (which
        # can't distinguish a genuine internal convergence from a
        # callback-triggered early stop: SPGBox returns the very same ierr==0
        # for both — spgbox_main.jl's `return SPGBoxResult(...,0,true)` on a
        # stall-detector hit).
        stall_reasons = String[]
        dmin_stalled_flag[] && push!(stall_reasons, "minimum distance plateaued")
        constraint_stalled_flag[] && push!(stall_reasons, "constraint violation plateaued")
        finish!(progress_meter)

        loop_end_reason = if tol_ok && const_ok
            "converged"
        elseif dmin_stalled_flag[] || constraint_stalled_flag[]
            "stalled (" * join(stall_reasons, ", ") * ")"
        elseif optresult.ierr == 0
            # Neither plateau flag fired, yet SPGBox still stopped with
            # ierr==0 short of tol_ok/const_ok: this can only be SPGBox's own
            # internal convergence criterion (small projected gradient) —
            # not a windowed plateau — typically because re-optimizing from
            # the exact same starting point produces so few iterations that
            # the plateau detectors never get the several consecutive
            # readings they need to fire, even though nothing is actually
            # still improving.
            "optimizer converged internally (small gradient), short of " *
                join(filter(!isnothing, [!tol_ok ? "the distance tolerance" : nothing,
                    !const_ok ? "the constraints" : nothing]), " and ")
        elseif optresult.ierr == 2
            "chunk function-evaluation budget (nfevalmax) reached"
        else
            "chunk iteration budget (maxit) reached"
        end

        @printf("\n  Packing loop ended: %s\n", loop_end_reason)
        @printf("  Function value from last loop: f = %10.5e\n", f_true_loop_end)
        @printf("  Best function value before: f = %10.5e\n", bestf_before_loop)
        @printf("  Improvement from best function value: %8.2f %%\n", fimprov)
        @printf("  Improvement within this loop: %8.2f %%\n", fimp_within_loop)
        @printf("  Minimum distance: %12.6f\n", dmin)
        @printf("  Maximum violation of the constraints: %10.5e\n", max_const)

        if tol_ok && const_ok
            println()
            println(hash_line)
            @printf("\n%s Success! \n", " "^32)
            @printf("%s Final objective function value: %10.5e\n", " "^13, f_true_loop_end)
            @printf("%s Minimum distance: %10.6f\n", " "^13, dmin)
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

        # Every atom's working radius always follows the same loop-indexed
        # radscale schedule toward its own floor (atom_radii_floor), every
        # loop, regardless of stall state — this is unconditional background
        # behavior, not a stall response. A molecule relocated by movebad!
        # below is not given its own separate, restarted schedule — its
        # atoms decay along with everyone else's.
        current_radscale = _radscale_at(loop + 1, packmol_system.radscale, radscale_decay_loops)
        for iat in eachindex(atom_radii)
            atom_radii[iat] = current_radscale * atom_radii_floor[iat]
        end
        # Whether to move bad molecules is decided purely at this outer-loop
        # granularity, on the true (radscale == 1.0) objective: if this
        # loop's own optimization improved `f_true_loop_end` over
        # `f_true_loop_start` by less than `movebad_tolerance` (10% by
        # default), that's this loop's own real progress falling short,
        # regardless of what atom_radii happened to be inflated to while
        # optimizing (movebad! itself still uses the current, possibly still
        # inflated atom_radii — via cl_system.fg.fmol, restored to that basis
        # just above — to pick and relocate molecules, and its own
        # exponential-probability-on-worst-fmol selection is unchanged).
        # The intra-chunk stall detectors above (`dmin_stalled_flag`/
        # `constraint_stalled_flag`) are no longer part of this decision —
        # they still exist purely to cut a chunk's SPGBox run short once it's
        # internally plateaued, so time isn't wasted grinding through the
        # rest of its maxit/nfevalmax budget.
        do_movebad = fimp_within_loop < T(100) * movebad_tolerance
        if do_movebad
            cm_min, cm_max = compute_cm_bounds(packmol_system)
            moved = movebad!(
                packmol_system, cl_system.fg.fmol, free_mol_indices, mol_structure_type, RNG;
                movefrac, precision, fmol_max_type,
                cm_lo_type=cm_min, cm_hi_type=cm_max,
                fixed_sys=movebad_fixed_sys, fixed_lo=movebad_fixed_lo, fixed_hi=movebad_fixed_hi,
                overlap_tol=movebad_overlap_tol,
                fg_output=movebad_fg_output, atom_positions=buffers.atom_positions, mol_iat_first,
            )
            if !isempty(moved)
                println("  Moved $(length(moved)) bad molecules randomly to new positions.")
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

    # For PBC: wrap each molecule's CM into the unit cell centered at
    # unitcell_center, and carry every atom of that molecule along by the
    # same offset (rigidly, via its CM) rather than wrapping each atom's
    # absolute position independently. Wrapping atoms independently lets a
    # molecule straddling a periodic boundary be torn in two (part of it
    # wrapped to the opposite face while the rest stays put), which can land
    # those wrapped atoms on top of whatever else sits there — visible as
    # overlapping atoms in the output even though the packing itself
    # converged. `write_output` recomputes atom positions from
    # `molecule_positions` below, so updating the CM here is what actually
    # takes effect. See the analogous comment on `_constraint_fg!` in
    # interatomic_distance_fg.jl for the same fix applied during optimization.
    if has_pbc
        center = packmol_system.unitcell_center
        for imol in eachindex(packmol_system.molecule_positions)
            mp = packmol_system.molecule_positions[imol]
            packmol_system.molecule_positions[imol] = MoleculePosition(
                wrap_to_center(mp.cm, packmol_system.unitcell, center), mp.angles
            )
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
    f_stall_tolerance = 0.05,
    f_true_loop_start::Union{Nothing,Real} = nothing,
    dmin_stalled_flag::Union{Nothing,Ref{Bool}} = nothing,
    constraint_stalled_flag::Union{Nothing,Ref{Bool}} = nothing,
    f_chunk_start::Union{Nothing,Ref} = nothing,
    f_true_start_progress_used::Union{Nothing,Ref{Bool}} = nothing,
    nfeval::Union{Nothing,Integer} = nothing,
    nfevalmax::Union{Nothing,Integer} = nothing,
    gnorm::Union{Nothing,Real} = nothing,
)
    dmin = min(cl_system.fg.dmin, cl_system.cutoff)
    max_const = cl_system.fg.max_constraint_penalty
    # nfeval/nfevalmax and gnorm are diagnostic only (not used in any
    # convergence/stall decision below): the progress bar and the
    # f/dmin/max_const values above only update once SPGBox *accepts* an
    # outer iteration, so a chunk stuck deep in one slow or failing internal
    # line search can look completely frozen for a long time even though
    # it's still working — nfeval climbing (against the nfevalmax budget
    # this whole chunk is capped at) shows that's what's happening, while a
    # blown-up or NaN gnorm points instead to a genuine numerical problem
    # (e.g. a degenerate/coincident atom pair) rather than just a hard,
    # slow-converging landscape.
    next!(progress_meter; showvalues = [
        (" Function value", f),
        (" Minimum distance", dmin),
        (" Maximum constraint violation", max_const),
        (" Function evaluations", isnothing(nfeval) ? nfeval : "$nfeval / $nfevalmax"),
        (" Projected gradient norm", gnorm),
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
    # (< f_stall_tolerance, e.g. 5%) over its own, much longer window
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
    # vetoed by whole-chunk progress: `f_chunk_start` starts out as the first
    # f this callback ever sees in this chunk, and as long as f has dropped
    # more than `f_stall_tolerance` from that reference, the chunk keeps
    # running regardless of what the trailing window shows — but the
    # reference is then advanced to the current f (see the re-baseline below)
    # every time that grants a reprieve, rather than staying pinned to the
    # chunk's opening value forever. Without that re-baseline this veto is a
    # one-way latch: a chunk that drops >5% in its first few dozen
    # iterations and then genuinely flatlines for the rest of its multi-
    # hundred-iteration budget would have `chunk_progress_significant` stuck
    # true for the remainder (the ratio against a fixed, already-cleared
    # opening value can only grow as f falls further, never revert), so the
    # dmin/const stall detectors below — and thus movebad!'s only trigger —
    # would never fire again for the rest of the chunk, no matter how long it
    # sits doing nothing; observed in practice as a chunk grinding silently
    # through its entire maxit/nfevalmax budget on a large system instead of
    # being cut short to let the outer loop relocate whatever molecule is
    # actually stuck. Re-baselining requires *renewed* progress to keep
    # earning the reprieve, which is what "genuine ongoing progress vetoes a
    # trailing-window plateau" was actually supposed to mean.
    if !isnothing(f_chunk_start) && !isnothing(f) && f_chunk_start[] == typemax(f_chunk_start[])
        f_chunk_start[] = oftype(dmin, f)
    end
    local_progress_significant = !isnothing(f_chunk_start) && !isnothing(f) &&
        f_chunk_start[] < typemax(f_chunk_start[]) && f_chunk_start[] > zero(f_chunk_start[]) &&
        (f_chunk_start[] - oftype(dmin, f)) / f_chunk_start[] > oftype(dmin, f_stall_tolerance)
    # In addition to the chunk-local check above (this chunk's own working,
    # atom_radii-inflated f, start vs now), also treat this chunk as making
    # significant progress, *once*, if it has already cleared
    # `f_stall_tolerance` relative to `f_true_loop_start` — the *true*
    # (floor-radii) objective at the exact position this loop started from,
    # fixed for the duration of this chunk and passed in by the caller. This
    # mirrors the outer loop's own "Improvement within this loop" check
    # (`fimp_within_loop`, computed once after the chunk ends, from
    # `f_true_loop_start` and the true floor-radii value at the chunk's end)
    # using the same reference — so a chunk the outer loop will end up
    # judging as having made real headway since it started isn't cut short
    # in here first purely for looking flat in working-scale terms, in the
    # first few iterations before `f_chunk_start`'s own window has enough
    # history to reflect that. `f` here is still this chunk's own working
    # objective, not the true one `f_true_loop_start` is measured in, so this
    # comparison mixes scales — `f` is not guaranteed to sit on either side
    # of `f_true_loop_start` in general; this is a useful, cheap proxy for
    # *early*-chunk progress, not an exact bound. Crucially, unlike
    # `local_progress_significant` above, this comparison is against a fixed
    # reference that only improves as f falls, so once true it would stay
    # true for the rest of the chunk no matter how long f then sits frozen —
    # exactly the one-way-latch failure mode described above, just with
    # `f_true_loop_start` playing the role of the never-advancing reference
    # instead of `f_chunk_start`. `f_true_start_progress_used` is the
    # re-baseline for *this* signal: once it has excused a plateau one time,
    # it's permanently spent for the rest of the chunk, so any further
    # "still improving" verdict must come from `local_progress_significant`
    # instead, which does require renewed progress each time.
    progress_vs_true_start = !isnothing(f_true_start_progress_used) && !f_true_start_progress_used[] &&
        !isnothing(f_true_loop_start) && !isnothing(f) &&
        f_true_loop_start > zero(f_true_loop_start) &&
        (f_true_loop_start - oftype(dmin, f)) / f_true_loop_start > oftype(dmin, f_stall_tolerance)
    chunk_progress_significant = local_progress_significant || progress_vs_true_start
    if chunk_progress_significant
        f_chunk_start[] = oftype(dmin, f)
        progress_vs_true_start && (f_true_start_progress_used[] = true)
    end
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
