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

#
# CellListMap cutoff for the interatomic-distance term. Every pair that can be
# penalized (d < r_i + r_j, under either the working radii `atom_radii` or the
# `atom_radii_floor` used by the "true objective" evaluations) must be within
# the cutoff. It is also at least `tolerance`, such that `dmin` (which is only
# known up to the cutoff, see `_format_dmin`) can always be compared with the
# tolerance in the convergence test. A larger cutoff only adds pairs that are
# computed and discarded (the cost of the pairwise loop grows with the number
# of pairs within the cutoff, not with the number of cells), so it is only
# grown to keep the number of cells below `_MAX_CELLS_PER_ATOM * natoms`.
#
const _MAX_CELLS_PER_ATOM = 2
function _distance_cutoff(atom_radii::Vector{T}, atom_radii_floor::Vector{T}, tolerance::T, volume::T, D::Int) where {T}
    rmax = max(maximum(atom_radii), maximum(atom_radii_floor))
    return _capped_cutoff(volume, max(2 * rmax, tolerance), _MAX_CELLS_PER_ATOM * length(atom_radii), D)
end

# Pairs farther than the cutoff are not computed, so a `dmin` equal to the
# cutoff only means that no pair is closer than it.
_format_dmin(dmin, cutoff) = dmin >= cutoff ? "> " * @sprintf("%.6f", cutoff) : @sprintf("%.6f", dmin)

#
# Unit cell (sides of an orthorhombic box) of the periodic system used by CellListMap
# for non-periodic systems: the region where the atoms are, inflated by 20%. The
# center-of-mass variables are not bounded in non-periodic systems, so molecules may
# move beyond this region during the optimization; the 20% slack absorbs that (a
# molecule would have to move about 10% of the box beyond the region before its atoms
# interact with the periodic images of the atoms on the opposite side).
#
function _nonperiodic_unitcell(packmol_system::PackmolSystem{D,T}, atom_positions::Vector{SVector{D,T}}) where {D,T}
    # Range of the current centers of the free molecules (`compute_cm_bounds` without
    # a `precision`, so not filtered by constraint satisfaction: no molecule is
    # projected onto its constraint region; not extrapolated: the actual range).
    cm_lo, cm_hi = compute_cm_bounds(packmol_system; extrapolate=false)
    # Overall bounds across all structure types
    all_lo = reduce((a,b) -> min.(a,b), cm_lo)
    all_hi = reduce((a,b) -> max.(a,b), cm_hi)
    # Margin: the greatest distance of an atom to the center of its molecule (the
    # reference coordinates of free structure types are centered at their geometric
    # center, the point the `cm` variables position), over the free structure types,
    # plus the tolerance. Fixed molecules are accounted for by the positions of their
    # atoms: a large fixed molecule (e.g. a solvated capsid) would otherwise inflate
    # the margin, and the cell volume, by its own radius.
    max_extent = zero(T)
    for st in packmol_system.structure_types
        st.fixed.fixed && continue
        for r in st.reference_coordinates
            max_extent = max(max_extent, norm(r))
        end
    end
    margin = max_extent + packmol_system.tolerance
    lo = all_lo .- margin
    hi = all_hi .+ margin
    iat = 0
    for st in packmol_system.structure_types
        nat = st.number_of_molecules * st.natoms
        if st.fixed.fixed
            for x in @view(atom_positions[iat+1:iat+nat])
                lo = min.(lo, x .- packmol_system.tolerance)
                hi = max.(hi, x .+ packmol_system.tolerance)
            end
        end
        iat += nat
    end
    box_size = hi - lo
    # Guard: if bounds are degenerate, fall back to bounding box of atom positions
    if any(box_size .≤ zero(T))
        lo, hi = compute_bounding_box(atom_positions)
        box_size = hi - lo
    end
    return T(1.2) * box_size
end

"""
    packmol(input_file::String; kargs...)

Read a Packmol input file, run the packing optimization, and write the output.
Returns the built `PackmolSystem` (see `packmol(::PackmolSystem)` for details).

Paths in the input file (structure files, `output`, `restart_from`)
are relative to the directory of the input file.

With no input file, a file explorer is opened to choose it (see [`Packmol.pick_file`](@ref)).
To disable the file explorer, set the environment variable `PACKMOL_GUI="false"`.
"""
function packmol(input_file::String; D::Int=3, T::DataType=Float64, kargs...)
    packmol_system = read_packmol_input(input_file; D, T)
    return packmol(packmol_system; kargs...)
end
"""
    packmol(; kargs...)

Choose a Packmol input file in a file browser, and run it.

The same as `packmol(input_file; kargs...)`, with the file chosen in the file browser of
[`Packmol.pick_file`](@ref) (disabled by setting the environment variable
`PACKMOL_GUI="false"`, in which case the file is asked for in the terminal).
"""
packmol(; kargs...) = packmol(pick_file(); kargs...)

"""
    packmol(packmol_system::PackmolSystem; kargs...)

Run the packing optimization on a `PackmolSystem`, mutating it in place (final
`molecule_positions`, ...) and returning that same object — not a bare `Bool` —
so callers can keep using it afterwards (e.g. to postprocess or re-`write_output`
the packed coordinates, or to convert their representation, as with
[`triclinic_to_dodecahedral`](@ref)).

The packing outcome is recorded in `packmol_system.status`: `:packing_ready` if
it converged within `nloop` loops (both `tolerance_precision` and
`constraint_precision` satisfied), `:failed_packing` if it did not, or
`:failed_initial_approximation` if the initial approximation could not satisfy the
constraints (written to the output file; the packing continues only if
`packmol_system.continue_failed_initial_packing` is set), or
`:not_packed` if `packmol_system.check_initial_placement` was set (initial approximation written,
no packing attempted).

**Global options.** Every global keyword of the input files (see [Global keywords](@ref)) can
be given as a keyword argument, with the name of the keyword: `tolerance`, `output`,
`filetype`, `seed`, `nloop`, `maxit`, `movefrac`, `radscale`, `radscale_decay_loops`,
`writeout`, `tolerance_precision`, `constraint_precision`, `nthreads`, `avoid_overlap`,
`connect`, `add_box_sides`, `add_amber_ter`, `amber_ter_preserve`, `periodic_boundary_style`,
`max_random_init`, `adjust_constraints_on_init`, `check_initial_placement`,
`continue_failed_initial_packing`, `restart_from`, `pbc` and `unitcell`. They are set in
`packmol_system`, replacing the values of the input file (or of the constructor), before
the packing. For example, `packmol(system; tolerance = 2.5, maxit = 50, output = "out.pdb")`.
The values are those of the input keywords: numbers, strings or symbols for the options
(`periodic_boundary_style = :octahedral`), `true`/`false` for the flags; `pbc` takes the
values of the input keyword (`pbc = [40, 40, 40]`, `pbc = (:octahedral, 0, 0, 0, 80)`), and
`unitcell` the six cell parameters or a unit cell matrix. A new `tolerance` sets the radius
of the atoms that have the default one (`tolerance / 2`).

The defaults: `nloop` `500`, `maxit` (optimizer iterations per loop) `30`, `movefrac`
(maximum fraction of the molecules of each type moved per loop) `0.05`, `writeout`
(minimum time, in seconds, between two writes of the output during the packing: the best
solution so far is written after a loop that improves it, if the end of the last write,
or the start of the packing, was at least `writeout` seconds before; the final solution is always
written) `60`, `nthreads` (maximum number of threads used) all the threads of the Julia
session, and `radscale_decay_loops` `20`: the number of loops in which the radii of the
atoms, multiplied at the start by the `radscale` of their structure type (see
[`PackmolSystem`](@ref); by default 5.0 for the free molecules of more than 20 atoms, 1.0
for the others), decay exponentially to their values (to within 0.1% of the initial
excess).

**Options of the run.**

`active_set` (default `true`): after the first loop, optimize only the
molecules with overlaps or constraint violations that may block the convergence, and the
molecules with an atom within `active_radius` (default `6 * tolerance`) of them, with the
other molecules frozen. Used in a loop when the atoms of the active molecules are at most
`active_max_fraction` (default `0.7`) of the atoms of the free molecules; otherwise the
whole system is optimized. Each loop reports the fractions of active molecules and atoms.
Short loops with the active subsystem pack large systems fastest: about 2.5x faster
than whole-system loops of 800 iterations for a nanoliposome of 600,000 atoms.

`rotation_scaling` (default `true`): the optimizer works on rotation variables scaled by
the inertia of each molecule (a unit change moves its atoms by 1 Å RMS), instead of its
Euler angles. With `false`, the Euler angles themselves, much slower for systems of
elongated molecules (lipids, polymers); for comparisons only.
"""
function packmol(packmol_system::PackmolSystem;
    output::Union{Nothing,AbstractString} = nothing,
    filetype::Union{Nothing,AbstractString,Symbol} = nothing,
    tolerance::Union{Nothing,Real} = nothing,
    radscale::Union{Nothing,Real} = nothing,
    discale::Union{Nothing,Real} = nothing,
    radscale_decay_loops::Union{Nothing,Integer} = nothing,
    writeout::Union{Nothing,Real} = nothing,
    tolerance_precision::Union{Nothing,Real} = nothing,
    constraint_precision::Union{Nothing,Real} = nothing,
    maxit::Union{Nothing,Integer} = nothing,
    nloop::Union{Nothing,Integer} = nothing,
    movefrac::Union{Nothing,Real} = nothing,
    nthreads::Union{Nothing,Integer} = nothing,
    seed::Union{Nothing,Integer} = nothing,
    max_random_init::Union{Nothing,Integer} = nothing,
    adjust_constraints_on_init::Union{Nothing,Bool} = nothing,
    add_amber_ter::Union{Nothing,Bool} = nothing,
    amber_ter_preserve::Union{Nothing,Bool} = nothing,
    add_box_sides::Union{Nothing,Bool} = nothing,
    periodic_boundary_style::Union{Nothing,AbstractString,Symbol} = nothing,
    connect::Union{Nothing,Bool} = nothing,
    avoid_overlap::Union{Nothing,Bool} = nothing,
    check_initial_placement::Union{Nothing,Bool} = nothing,
    continue_failed_initial_packing::Union{Nothing,Bool} = nothing,
    restart_from::Union{Nothing,AbstractString} = nothing,
    pbc = nothing,
    unitcell = nothing,
    kargs...,
)
    # The global options of the input files, set in the system
    _set_global_options!(packmol_system;
        output, filetype, tolerance, radscale, discale, radscale_decay_loops, writeout,
        tolerance_precision, constraint_precision, maxit, nloop, movefrac, nthreads, seed,
        max_random_init, adjust_constraints_on_init, add_amber_ter, amber_ter_preserve,
        add_box_sides, periodic_boundary_style, connect, avoid_overlap, check_initial_placement,
        continue_failed_initial_packing, restart_from, pbc, unitcell,
    )
    ntasks = _ntasks[]
    _ntasks[] = something(packmol_system.nthreads, 0)
    try
        return _packmol(packmol_system; kargs...)
    finally
        _ntasks[] = ntasks
    end
end

# The objective function of the scaled variables `y` (see rigid_body/rotation_metric.jl):
# evaluated at the Euler angles (in `x_unscaled`), with the gradient relative to the
# scaled angles.
function _scaled_fg!(g, y, x_unscaled, rot_scaling, cl_system, packmol_system::PackmolSystem{D,T}, atom_positions, free_mol_indices, atom_radii) where {D,T}
    _unscale_angles!(reinterpret(MoleculePosition{D,T}, x_unscaled), reinterpret(MoleculePosition{D,T}, y), rot_scaling)
    f = fg!(g, x_unscaled, cl_system, packmol_system, atom_positions, free_mol_indices, atom_radii)
    _scale_angle_gradient!(reinterpret(MoleculePosition{D,T}, g), rot_scaling)
    return f
end

# The scaling of the rotations exists for three-dimensional systems only
_rotation_scaling(packmol_system::PackmolSystem{3}) = RotationScaling(packmol_system)
_rotation_scaling(::PackmolSystem) = nothing

# The local solver of the packing loop: SPGBox. A separate function so that the solver
# experiments of `dev/solvers/` can replace it by redefining this method. The bounds are
# positional: they are `nothing` without PBC, and a keyword call with arguments of a
# `Union` type is not resolved by `juliac --trim`.
_local_solver(fg!, x, lower, upper; callback, vaux, nitmax, nfevalmax) =
    spgbox!(fg!, x; callback, vaux, nitmax, nfevalmax, lower, upper)

function _packmol(
    packmol_system::PackmolSystem{D,T};
    parallel::Bool=true,
    iprint::Int=10,
    nloop::Union{Nothing,Int}=nothing,
    maxit::Union{Nothing,Int}=nothing,
    movefrac::Union{Nothing,Real}=nothing,
    radscale_decay_loops::Union{Nothing,Int}=nothing,
    writeout::Union{Nothing,Real}=nothing,
    movebad_tolerance::T=T(0.10),
    seed::Int=packmol_system.seed,
    restart::Bool=false,
    rotation_scaling::Bool=true,
    active_set::Bool=true,
    active_radius::Real=6 * packmol_system.tolerance,
    active_max_fraction::Real=0.7,
) where {D,T}
    # Keyword arguments take precedence over the input file's `nloop`/`maxit`
    nloop = something(nloop, packmol_system.nloop, 500)
    maxit = something(maxit, packmol_system.max_iter, 30)
    movefrac = T(something(movefrac, packmol_system.movefrac, 0.05))
    # Initialize RNG and molecule positions. Negative seed (Fortran Packmol's
    # `seed -1`, used by the Recipes as their default) means "pick a random
    # seed" rather than a literal RNG seed — passing it straight to
    # Xoshiro(seed) throws a DomainError on Julia 1.10 (LTS), since negative
    # integers there hit make_seed's "n must be non-negative" path.
    RNG = seed < 0 ? Random.Xoshiro() : Random.Xoshiro(seed)

    tstart = time()    

    # Print title
    _println()
    _println(hash_line)
    _println("  PACKMOL - Packing optimization for the automated generation of")
    _println("  starting configurations for molecular dynamics simulations.")
    _println()
    @_printf("%62s\n", "Version $_packmol_version ")
    _println(hash_line)
    _println()
    _check_fixed_inside_cell(packmol_system)
    _check_restart_files(packmol_system)

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
    # A *whole-system* restart source is whatever a prior
    # `write_output` call wrote for the entire
    # system, so it lists every molecule of every structure type, fixed
    # ones included, interleaved in structure-declaration order — not just
    # the free ones. `segments` below therefore spans every structure type
    # (not just the non-fixed ones), with `extract` (the tuple's 4th field)
    # marking which ones to actually read a position back from; a fixed
    # segment's atoms still have to be counted and skipped over so the
    # expected total and the file's actual layout line up (see
    # `_restart_positions_from_atoms`/`_restart_positions`), even though its
    # own position is discarded.
    _center_reference_coordinates!(packmol_system)
    restart_from = packmol_system.restart_from
    if !isnothing(restart_from)
        _println("  Restarting all free molecules from: ", _restart_source_name(restart_from))
        segments = Tuple{Int,Int,Vector{SVector{D,T}},Bool}[
            (_ncopies(st), st.natoms, st.reference_coordinates, !st.fixed.fixed)
            for st in packmol_system.structure_types
        ]
        positions = _restart_positions(restart_from, segments, MoleculePosition{D,T})
        for (k, imol) in enumerate(free_mol_indices)
            packmol_system.molecule_positions[imol] = positions[k]
        end
        restart = true
    end

    # Pre-optimization: set initial approximation (random placement + constraint fitting)
    initial_converged = set_initial_approximation!(packmol_system, free_mol_indices, RNG; restart, nloop, buffers)

    imol_offset = 0
    for st in packmol_system.structure_types
        st_restart_from = st.restart_from
        if !st.fixed.fixed && !isnothing(st_restart_from)
            _println("  Restarting structure type '", basename(st.filename), "' from: ", _restart_source_name(st_restart_from))
            # A per-structure-type restart source only ever contains this
            # one type's own molecules (e.g. the coordinates of that type from
            # an earlier run), so its segment list is just itself — always
            # extracted, nothing to skip.
            segments = [(st.number_of_molecules, st.natoms, st.reference_coordinates, true)]
            positions = _restart_positions(st_restart_from, segments, MoleculePosition{D,T})
            packmol_system.molecule_positions[imol_offset+1:imol_offset+st.number_of_molecules] .= positions
        end
        imol_offset += st.number_of_molecules
    end

    # The constraints could not be satisfied even without the distances between the
    # molecules: they are probably wrong (incompatible, or too tight for the molecules),
    # and packing would not converge. Stop, unless asked to continue.
    if !initial_converged && !packmol_system.check_initial_placement && !packmol_system.continue_failed_initial_packing
        _println()
        _println("  ERROR: the initial approximation does not satisfy the constraints (or overlaps")
        _println("  the fixed molecules). The constraints are probably wrong: incompatible with each")
        _println("  other, or too tight for the molecules. The initial approximation is written to the")
        _println("  output file, for inspection. To pack anyway, use the keyword")
        _println("  `continue_failed_initial_packing`.")
        if !isempty(packmol_system.output_file)
            write_output(packmol_system)
        end
        packmol_system.status = :failed_initial_approximation
        return packmol_system
    end

    # check mode: write the initial approximation and return
    if packmol_system.check_initial_placement
        _println("  Check mode: writing initial approximation to $(packmol_system.output_file)")
        if !isempty(packmol_system.output_file)
            write_output(packmol_system)
        end
        packmol_system.status = :not_packed
        return packmol_system
    end

    # Compute initial atom positions (reuse pre-allocated buffer)
    natoms = length(packmol_system.atoms)
    atom_positions = buffers.atom_positions
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)

    # Determine unit cell for CellListMap: a static matrix in both cases (diagonal
    # without PBC), so that the type of `cl_system` does not depend on `has_pbc`.
    has_pbc = !isnothing(packmol_system.unitcell)
    if has_pbc
        # PBC mode: use the actual unit cell
        unitcell = SMatrix{D,D,T}(packmol_system.unitcell)
    else
        unitcell = SMatrix{D,D,T}(Diagonal(_nonperiodic_unitcell(packmol_system, atom_positions)))
    end

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
    # The scale of the radii of each atom decays from that of its structure type
    # (`_radscale`, at loop == 0) to 1.0 by loop == radscale_decay_loops (default: 20).
    decay_loops = T(something(radscale_decay_loops, packmol_system.radscale_decay_loops, 20))
    atom_radscale = T[_radscale(packmol_system, packmol_system.structure_types[a.structure_type_index])
                      for a in packmol_system.atoms]
    # Minimum time between two writes of the output file during the packing
    write_interval = Float64(something(writeout, packmol_system.writeout, 60))
    last_write = time()
    atom_radii = atom_radscale .* atom_radii_floor

    # CellListMap cutoff: set from the radii in use (see `_distance_cutoff`),
    # and updated whenever `atom_radii` shrinks along the radscale schedule.
    tol = packmol_system.tolerance
    volume = abs(det(unitcell))
    cutoff = _distance_cutoff(atom_radii, atom_radii_floor, tol, volume, D)

    # Set up CellListMap
    fg_output = InteratomicDistanceFG{D,T}(packmol_system)
    _println("  Total number of atoms: ", natoms)
    _println("  Number of free molecules: ", nfree)
    _println("  Number of variables: ", nfree * 2 * D)
    _println("  Number of threads: ", _nbatches())
    cl_system = ParticleSystem(
        xpositions=atom_positions,
        unitcell=unitcell,
        cutoff=cutoff,
        output=fg_output,
        output_name=Val(:fg),
        parallel=parallel && _nbatches() > 1,
        nbatches=_cl_nbatches(),
        validate_coordinates=_validate_positions,
    )

    # Set up optimization variables: only free molecule DOFs (reuse pre-allocated buffer)
    _wrap_molecule_positions!(packmol_system, free_mol_indices)
    x = buffers.x
    x_mol = reinterpret(MoleculePosition{D,T}, x)
    for (k, imol) in enumerate(free_mol_indices)
        x_mol[k] = packmol_system.molecule_positions[imol]
    end
    auxvecs = buffers.vaux

    # Scaling of the rotations (see rigid_body/rotation_metric.jl): in each chunk, the
    # optimizer works on the angles a′ = Lᵀ·a of each molecule, with which a unit change
    # moves its atoms by 1 Å RMS. `x_unscaled` holds the variables in the original
    # (Euler angles) form, for the objective function. Three-dimensional packings only.
    scale_rotations = rotation_scaling && D == 3
    rot_scaling = _rotation_scaling(packmol_system)
    x_unscaled = similar(x)

    # Placement region for movebad! randomization
    mol_structure_type = _build_mol_structure_type(packmol_system)
    mol_iat_first = _build_mol_iat_first(packmol_system)
    precision = packmol_system.tolerance_precision
    # Worst fmol ever observed per structure type, across every movebad! call
    # over the whole run (see movebad.jl for why this must persist rather
    # than being recomputed fresh each call).
    fmol_max_type = zeros(T, length(packmol_system.structure_types))
    # The molecules that do not satisfy the convergence criterion (the bad molecules of
    # the active subsystem and of movebad!), at the start and at the end of each loop
    mol_unconverged = zeros(Bool, packmol_system.nmols)

    # Fixed-structure overlap-check system and constraint-only scratch for
    # movebad! (see movebad.jl): built once here, since the fixed atoms never
    # move over the course of the run, and reused on every relocation.
    movebad_overlap_tol = _global_radscale(packmol_system) * packmol_system.tolerance
    movebad_fixed_sys, movebad_fixed_lo, movebad_fixed_hi =
        _build_overlap_check_system(packmol_system, movebad_overlap_tol; avoid_overlap=true)
    movebad_fg_output = InteratomicDistanceFG{D,T}(packmol_system)

    # Under PBC, the centers of mass are hard-bounded — to exactly one
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
    # several periods and having to claw back. The bounds are closed, and the
    # upper faces belong to the cell: centers of mass inside the closed cell are
    # not wrapped when evaluating constraints (`_wrap_cm`). `nothing` without
    # PBC, so SPGBox skips the bound-checking overhead entirely.
    lower, upper = if has_pbc
        # Axis-aligned half-extent of the (possibly triclinic) periodic
        # cell — exact for an orthorhombic (diagonal) unitcell, and a
        # (safe, if not perfectly tight) enclosing box for a sheared one,
        # and for the rhombic dodecahedron or truncated octahedron of
        # dodecahedral and octahedral cells (see `_cm_half_extent`).
        half_extent = _cm_half_extent(packmol_system)
        cm_lo = packmol_system.unitcell_center .- half_extent
        cm_hi = packmol_system.unitcell_center .+ half_extent
        ang_lo = SVector{D,T}(ntuple(_ -> T(-Inf), D))
        ang_hi = SVector{D,T}(ntuple(_ -> T(Inf), D))
        reinterpret(T, fill(MoleculePosition(cm_lo, ang_lo), nfree)),
            reinterpret(T, fill(MoleculePosition(cm_hi, ang_hi), nfree))
    else
        nothing, nothing
    end

    # Atoms of the free molecules (for the report of the active subsystem)
    nfree_atoms = sum(packmol_system.structure_types[mol_structure_type[imol]].natoms for imol in free_mol_indices; init=0)

    # Outer packing loop (following Fortran Packmol gencanloop):
    # Each iteration runs a short optimization, evaluates per-molecule
    # contributions, and randomly re-places the worst molecules.
    _println()
    _println(dash_line)
    _println("  Packing $nfree free molecules ($(packmol_system.nmols) total)...")
    _println(dash_line)
    # Evaluate and print initial function value
    g0 = similar(x)
    # Objective function at the initial point, on the same unscaled
    # (radscale == 1.0) basis as every other cross-loop comparison below.
    f0_true = fg!(g0, x, cl_system, packmol_system, atom_positions, free_mol_indices, atom_radii_floor)
    @_printf("  Objective function at initial point: %10.5e\n", f0_true)
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

        _println()
        _println(dash_line)
        @_printf("  Starting packing loop: %8d\n", loop)
        atom_radii_lo, atom_radii_hi = extrema(atom_radii)
        @_printf("  Atom radii in this loop: [ %8.4f - %8.4f ]\n", atom_radii_lo, atom_radii_hi)
        _println()

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
        # (recording the molecules that do not satisfy the convergence criterion: the
        # bad molecules of the active subsystem)
        cl_system.fg.count_overlaps[] = true
        f_true_loop_start = fg!(g0, x, cl_system, packmol_system, atom_positions, free_mol_indices, atom_radii_floor)
        cl_system.fg.count_overlaps[] = false
        _mark_unconverged!(mol_unconverged, cl_system.fg)

        # Run a short optimization (maxit iterations per loop)
        progress_meter = ProgressBar(maxit; desc="  Iterations: ", barlen=46)
        # On an active subset of the system (see active_set.jl): the bad molecules and
        # their neighbors, from the contributions to the objective function at the
        # start of the loop (`cl_system.fg.fmol`, of the evaluation above)
        sub, nactive, nactive_atoms = if active_set && loop > 0
            _active_subsystem(packmol_system, mol_unconverged, free_mol_indices, cl_system.xpositions,
                unitcell, cl_system.cutoff, mol_structure_type, mol_iat_first;
                active_radius=T(active_radius), max_fraction=T(active_max_fraction))
        else
            nothing, 0, 0
        end
        if active_set && loop > 0
            @_printf("  Active molecules: %d of %d (%.2f%%)\n", nactive, nfree, 100 * nactive / nfree)
            @_printf("  Active atoms: %d of %d (%.2f%%)\n", nactive_atoms, nfree_atoms, 100 * nactive_atoms / nfree_atoms)
            _println(isnothing(sub) ? "  Above active_max_fraction: optimizing the whole system" : "  Optimizing the active subsystem")
            _println()
        end
        optresult = if isnothing(sub)
            if scale_rotations
                _update_rotation_scaling!(rot_scaling, reinterpret(MoleculePosition{D,T}, x))
                _scale_angles!(reinterpret(MoleculePosition{D,T}, x), rot_scaling)
            end
            fg_closure = (g, y) -> scale_rotations ?
                _scaled_fg!(g, y, x_unscaled, rot_scaling, cl_system, packmol_system, atom_positions, free_mol_indices, atom_radii) :
                fg!(g, y, cl_system, packmol_system, atom_positions, free_mol_indices, atom_radii)
            res = _local_solver(
                fg_closure,
                x, lower, upper;
                callback=(result) -> packmol_callback(cl_system, tol, iprint,
                    packmol_system.tolerance_precision, packmol_system.constraint_precision, progress_meter;
                    f=result.f, nfeval=result.nfeval, nfevalmax=10 * maxit, gnorm=result.gnorm,
                ),
                vaux=auxvecs,
                nitmax=maxit,
                nfevalmax=10 * maxit,
            )
            # Back to the Euler angles
            if scale_rotations
                _unscale_angles!(reinterpret(MoleculePosition{D,T}, x), reinterpret(MoleculePosition{D,T}, x), rot_scaling)
            end
            res
        else
            na = length(sub.mols)
            x_mol_full = reinterpret(MoleculePosition{D,T}, x)
            x_act = Vector{T}(undef, 2 * D * na)
            xa_mol = reinterpret(MoleculePosition{D,T}, x_act)
            for k in 1:na
                xa_mol[k] = x_mol_full[sub.kfree[k]]
            end
            lower_act, upper_act = if isnothing(lower) || isnothing(upper)
                nothing, nothing
            else
                collect(reinterpret(T, reinterpret(MoleculePosition{D,T}, lower)[sub.kfree])),
                    collect(reinterpret(T, reinterpret(MoleculePosition{D,T}, upper)[sub.kfree]))
            end
            rs_act = scale_rotations ? _subset(rot_scaling, sub.kfree) : nothing
            if !isnothing(rs_act)
                _update_rotation_scaling!(rs_act, xa_mol)
                _scale_angles!(xa_mol, rs_act)
            end
            xu_act = similar(x_act)
            fg_act = (g, y) -> if isnothing(rs_act)
                _active_fg!(g, y, sub, packmol_system, atom_radii)
            else
                _unscale_angles!(reinterpret(MoleculePosition{D,T}, xu_act), reinterpret(MoleculePosition{D,T}, y), rs_act)
                fa = _active_fg!(g, xu_act, sub, packmol_system, atom_radii)
                _scale_angle_gradient!(reinterpret(MoleculePosition{D,T}, g), rs_act)
                fa
            end
            res = _local_solver(
                fg_act,
                x_act, lower_act, upper_act;
                callback=(result) -> packmol_callback(sub.cl, tol, iprint,
                    packmol_system.tolerance_precision, packmol_system.constraint_precision, progress_meter;
                    f=result.f, nfeval=result.nfeval, nfevalmax=10 * maxit, gnorm=result.gnorm,
                ),
                vaux=SPGBox.VAux(x_act, zero(T)),
                nitmax=maxit,
                nfevalmax=10 * maxit,
            )
            if !isnothing(rs_act)
                _unscale_angles!(xa_mol, xa_mol, rs_act)
            end
            for k in 1:na
                x_mol_full[sub.kfree[k]] = xa_mol[k]
            end
            res
        end

        # Update molecule positions with optimized values
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(free_mol_indices)
            packmol_system.molecule_positions[imol] = x_mol[k]
        end
        _wrap_molecule_positions!(packmol_system, free_mol_indices)

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
        # (also counting the overlapping pairs, for the report of the loop, and recording
        # the molecules that do not satisfy the convergence criterion: the bad molecules
        # of movebad!)
        cl_system.fg.count_overlaps[] = true
        f_true_loop_end = fg!(g0, x, cl_system, packmol_system, atom_positions, free_mol_indices, atom_radii_floor)
        cl_system.fg.count_overlaps[] = false
        _mark_unconverged!(mol_unconverged, cl_system.fg)
        overlaps = copy(cl_system.fg.overlaps)
        dmin = min(cl_system.fg.dmin, cl_system.cutoff)
        max_overlap = cl_system.fg.max_overlap
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

        # Check convergence: both tolerance and constraint precisions must be satisfied.
        # The overlap of each pair is measured with the radii of its atoms (by default
        # tolerance/2, so that the criterion is `tolerance - dmin < tolerance_precision`).
        tol_ok = max_overlap < packmol_system.tolerance_precision
        const_ok = max_const < packmol_system.constraint_precision

        finish!(progress_meter)
        loop_end_reason = if tol_ok && const_ok
            "converged"
        elseif optresult.ierr == 0
            # SPGBox's own convergence criterion (small projected gradient)
            "small gradient, " *
                (!tol_ok && !const_ok ? "tolerance and constraints not satisfied" :
                    !tol_ok ? "tolerance not satisfied" : "constraints not satisfied")
        elseif optresult.ierr == 2
            "function-evaluation budget (10 maxit) reached"
        else
            "iteration budget (maxit) reached"
        end

        @_printf("\n  Packing loop ended: %s\n", loop_end_reason)
        _print_overlaps(packmol_system, overlaps) && _println()
        @_printf("  Function value from last loop: f = %10.5e\n", f_true_loop_end)
        @_printf("  Best function value before: f = %10.5e\n", bestf_before_loop)
        @_printf("  Improvement from best function value: %8.2f %%\n", fimprov)
        @_printf("  Improvement within this loop: %8.2f %%\n", fimp_within_loop)
        @_printf("  Minimum distance: %12s\n", _format_dmin(dmin, cl_system.cutoff))
        @_printf("  Maximum violation of the constraints: %10.5e\n", max_const)

        if tol_ok && const_ok
            _println()
            _println(hash_line)
            @_printf("\n%32s Success! \n", "")
            @_printf("%13s Final objective function value: %10.5e\n", "", f_true_loop_end)
            @_printf("%13s Minimum distance: %10s\n", "", _format_dmin(dmin, cl_system.cutoff))
            @_printf("%13s Maximum violation of the constraints: %10.5e\n", "", max_const)
            _println()
            _println(dash_line)
            _println()
            _println("              If Packmol contributed to your work, please cite: ")
            _println()
            _println("            L. Martinez, R. Andrade, E. G. Birgin, J. M. Martinez, ")
            _println("          PACKMOL: A package for building initial configurations for")
            _println("                   molecular dynamics simulations. ")
            _println("         Journal of Computational Chemistry, 30(13) pp. 2157-2164, 2009.")
            _println("                  https://doi.org/10.1002/jcc.21224")
            _println()
            _println(hash_line)
            converged = true
            break
        end

        _println()
        _println(dash_line)

        # Write best solution so far to output file (at most every `write_interval`
        # seconds after the end of the last write: writing a large system can take
        # longer than a loop)
        if improved && !isempty(packmol_system.output_file) && time() - last_write >= write_interval
            saved_positions = packmol_system.molecule_positions
            packmol_system.molecule_positions = best_positions
            write_output(packmol_system)
            packmol_system.molecule_positions = saved_positions
            last_write = time()
            _println()
            _println("  Current solution written to file: ", packmol_system.output_file)
        end

        # Every atom's working radius always follows the same loop-indexed
        # radscale schedule toward its own floor (atom_radii_floor), every
        # loop, regardless of stall state — this is unconditional background
        # behavior, not a stall response. A molecule relocated by movebad!
        # below is not given its own separate, restarted schedule — its
        # atoms decay along with everyone else's.
        for iat in eachindex(atom_radii)
            atom_radii[iat] = _radscale_at(loop + 1, atom_radscale[iat], decay_loops) * atom_radii_floor[iat]
        end
        CellListMap.update!(cl_system; cutoff=_distance_cutoff(atom_radii, atom_radii_floor, tol, volume, D))
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
        # Never after loop 0: it starts from the initial placement, and with inflated radii
        # (radscale) it can increase the true objective function while it spreads the
        # molecules apart, which is no sign of molecules that need to be moved.
        do_movebad = loop > 0 && fimp_within_loop < T(100) * movebad_tolerance
        if do_movebad
            cm_min, cm_max = compute_cm_bounds(packmol_system)
            moved = movebad!(
                packmol_system, cl_system.fg.fmol, free_mol_indices, mol_structure_type, RNG;
                movefrac, precision, fmol_max_type, unconverged=mol_unconverged,
                cm_lo_type=cm_min, cm_hi_type=cm_max,
                fixed_sys=movebad_fixed_sys, fixed_lo=movebad_fixed_lo, fixed_hi=movebad_fixed_hi,
                overlap_tol=movebad_overlap_tol,
                fg_output=movebad_fg_output, atom_positions=buffers.atom_positions, mol_iat_first,
            )
            if !isempty(moved)
                _println("  Moved $(length(moved)) bad molecules randomly to new positions.")
            end
        end
        # Re-pack optimizer variables from (possibly moved) molecule positions
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(free_mol_indices)
            x_mol[k] = packmol_system.molecule_positions[imol]
        end
    end

    if !converged
        _println()
        _println(hash_line)
        @_printf("  WARNING: packing did not converge after %d loops (best f = %.4e)\n", nloop, bestf)
        _println(hash_line)
    end
    packmol_system.status = converged ? :packing_ready : :failed_packing

    # Restore best molecule positions
    copyto!(packmol_system.molecule_positions, best_positions)

    # Note: `packmol_system.molecule_positions` is left exactly as the
    # optimizer produced it — under PBC, a molecule's center of mass may sit
    # outside the unit cell (an unwrapped, but physically equivalent, periodic
    # image), which is fine since every distance/constraint evaluation during
    # optimization was already PBC-aware. Wrapping into the cell (rigidly, by
    # CM, so a molecule is never torn in two across a boundary — see
    # `get_atoms`) happens on demand instead, in `get_atoms`/`write_output`,
    # not by mutating this array here.

    # Write output file if specified
    if !isempty(packmol_system.output_file)
        write_output(packmol_system)
        _println()
        _println("  Solution written to file: ", packmol_system.output_file)
        converged || _write_forced_output(packmol_system, free_mol_indices, RNG, buffers, nloop)
    end

    _println()
    _println(dash_line)
    tend = time()
    @_printf("  Running time: %12.4f seconds.\n", tend - tstart)
    _println(dash_line)
    _println()

    return packmol_system
end

#
# Aesthetic line constants (matching Fortran Packmol)
#
const dash_line = repeat('-', 80)
const hash_line = repeat('#', 80)
const _packmol_version = string(pkgversion(@__MODULE__))

#
# When packing fails, write (as Fortran Packmol) the "forced" solution to
# `<output>_FORCED`: the best solution, after an adjustment that satisfies the
# constraints regardless of overlaps. The best solution itself is kept in the
# system (and in the regular output and restart files).
#
function _write_forced_output(packmol_system::PackmolSystem, free_mol_indices::Vector{Int}, RNG, buffers::MemoryBuffers, nloop::Int)
    best_positions = copy(packmol_system.molecule_positions)
    _println()
    _println("  Forcing the solution to fit the constraints...")
    cm_min, cm_max = compute_cm_bounds(packmol_system)
    adjust_constraints!(packmol_system, free_mol_indices, RNG;
        nloop, cm_lo_type=cm_min, cm_hi_type=cm_max, buffers,
    )
    forced_file = packmol_system.output_file * "_FORCED"
    write_output(packmol_system; output_file=forced_file)
    _println("  The forced point was written to the output file: ", forced_file)
    copyto!(packmol_system.molecule_positions, best_positions)
    return forced_file
end

#
# SPGBox callback: progress bar, and stop on convergence
#
function packmol_callback(
    cl_system, tol, iprint, tolerance_precision, constraint_precision, progress_meter;
    f::Union{Nothing,Real} = nothing,
    nfeval::Union{Nothing,Integer} = nothing,
    nfevalmax::Union{Nothing,Integer} = nothing,
    gnorm::Union{Nothing,Real} = nothing,
)
    dmin = min(cl_system.fg.dmin, cl_system.cutoff)
    max_const = cl_system.fg.max_constraint_penalty
    # nfeval/nfevalmax and gnorm are diagnostic only: the progress bar and the
    # values above only update once SPGBox accepts an iteration, so nfeval shows
    # a slow line search, and a blown-up or NaN gnorm a numerical problem.
    next!(progress_meter; showvalues = (
        ("Function value", isnothing(f) ? "" : @sprintf("%.5e", f)),
        ("Minimum distance", _format_dmin(dmin, cl_system.cutoff)),
        ("Maximum constraint violation", @sprintf("%.5e", max_const)),
        ("Function evaluations", isnothing(nfeval) ? "" : "$nfeval / $nfevalmax"),
        ("Projected gradient norm", isnothing(gnorm) ? "" : @sprintf("%.5e", gnorm)),
    ))
    tol_ok = cl_system.fg.max_overlap < tolerance_precision
    const_ok = max_const < constraint_precision
    return tol_ok && const_ok
end

#
# Report of the pairs of atoms overlapping by more than tolerance_precision (the
# convergence criterion), for the pairs of structure types with the most of them
#
function _print_overlaps(packmol_system::PackmolSystem, overlaps::Matrix{Int}; nshow::Int=6)
    total = sum(overlaps)
    total == 0 && return false
    @_printf("  Overlapping pairs of atoms: %d, between molecules of types:\n", total)
    pairs = Tuple{Int,Int,Int}[]
    for j in axes(overlaps, 2), i in 1:j
        overlaps[i, j] > 0 && push!(pairs, (overlaps[i, j], i, j))
    end
    sort!(pairs; rev=true)
    for k in 1:min(nshow, length(pairs))
        n, i, j = pairs[k]
        @_printf("    %10d  %s - %s\n", n, _type_label(packmol_system, i), _type_label(packmol_system, j))
    end
    if length(pairs) > nshow
        @_printf("    %10d  (other pairs of types)\n", sum(first(p) for p in pairs[nshow+1:end]))
    end
    return true
end
function _type_label(packmol_system::PackmolSystem, itype::Int)
    st = packmol_system.structure_types[itype]
    label = string(itype, ": ", basename(st.filename))
    return st.fixed.fixed ? label * " (fixed)" : label
end
