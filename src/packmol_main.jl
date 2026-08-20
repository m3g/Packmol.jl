"""
    packmol(input_file::String; kargs...)

Read a Packmol input file, run the packing optimization, and write the output.
"""
function packmol(input_file::String; D::Int=3, T::DataType=Float64, kargs...)
    packmol_system = read_packmol_input(input_file; D, T)
    packmol(packmol_system; kargs...)
    return nothing
end


"""
    packmol(packmol_system::PackmolSystem; kargs...)

Run the packing optimization on a `PackmolSystem`.

"""
function packmol(
    packmol_system::PackmolSystem{D,T};
    parallel::Bool=true,
    iprint::Int=10,
    nloop::Int=200,
    maxit::Union{Nothing,Int}=nothing,
    movefrac::T=T(0.05),
    seed::Int=packmol_system.seed,
    restart::Bool=false,
    optimizer::Symbol=:spgbox,
) where {D,T}
    optimizer in (:spgbox, :gencan) ||
        error("unknown optimizer :$optimizer (expected :spgbox or :gencan)")
    if optimizer === :gencan && T !== Float64
        error("optimizer=:gencan requires T === Float64 (GENCAN is hardcoded double precision), got T = $T")
    end
    # GENCAN's default per-chunk iteration budget matches the original
    # Fortran Packmol default (getinp.f90: maxit = 20); SPGBox keeps its
    # existing default of 200. Either can be overridden explicitly.
    maxit = something(maxit, optimizer === :gencan ? 20 : 200)
    # Initialize RNG and molecule positions
    RNG = Random.Xoshiro(seed)

    tstart = time()    

    # Print title
    _version = pkgversion(Packmol)
    println()
    println(hash_line)
    println(" PACKMOL - Packing optimization for the automated generation of")
    println(" starting configurations for molecular dynamics simulations.")
    println(" ")
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

    # Pre-optimization: set initial approximation (random placement + constraint fitting)
    set_initial_approximation!(packmol_system, free_mol_indices, RNG; restart, buffers)

    # check mode: write the initial approximation and return
    if packmol_system.check
        println("Check mode: writing initial approximation to $(packmol_system.output_file)")
        if !isempty(packmol_system.output_file)
            write_output(packmol_system)
        end
        return nothing
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
    println(" Total number of atoms: ", natoms)
    println(" Number of free molecules: ", nfree)
    println(" Number of variables: ", nfree * 2 * D)
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

    # Outer packing loop (following Fortran Packmol gencanloop):
    # Each iteration runs a short optimization, evaluates per-molecule
    # contributions, and randomly re-places the worst molecules.
    println()
    println(dash_line)
    println(" Packing $nfree free molecules ($(packmol_system.nmols) total)...")
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
    @printf(" Objective function at initial point: %10.5e\n", f0)
    bestf = typemax(T)
    flast = typemax(T)
    converged = false
    best_positions = copy(packmol_system.molecule_positions)
    for loop in 0:nloop

        println()
        println(dash_line)
        @printf(" Starting packing loop: %8d\n", loop)
        @printf(" Tolerance in this loop: %8.4f\n", radscale * tol)
        println()

        # Run a short optimization (maxit iterations per loop)
        progress_meter = Progress(maxit; desc=" Iterations: ", barlen=47)
        fg_closure = (g, x) -> fg!(g, x, cl_system, packmol_system, atom_positions, free_mol_indices, radscale)
        # Fresh stall detector per chunk: plateau state from the *previous*
        # chunk isn't meaningful here, since movebad! (or the initial
        # placement) just changed the starting point.
        stall_detector = StallDetector{T}()
        optresult = if optimizer === :spgbox
            spgbox!(
                fg_closure,
                x;
                callback=(result) -> packmol_callback(result, cl_system, tol, iprint,
                    packmol_system.tolerance_precision, packmol_system.constraint_precision, progress_meter;
                    stall_detector, stall_tolerance=packmol_system.stall_tolerance,
                ),
                vaux=auxvecs,
                nitmax=maxit,
                nfevalmax=10 * maxit,
            )
        else # optimizer === :gencan
            # GENCAN's callback (unlike SPGBox's) can't signal early-exit —
            # it runs its full maxit-iteration budget (or until its own
            # internal epsgpsn criterion triggers). This is fine: the outer
            # loop below still re-checks Packmol's own tol_ok/const_ok
            # convergence after every chunk regardless of which optimizer
            # produced it. `packmol_callback`'s first argument (`spgresult`)
            # is unused in its body, so it's safe to reuse here with a
            # dummy value just for the progress-bar/showvalues side effect.
            # Stall detection is skipped too (stall_detector=nothing): it
            # relies on `spgresult.f` at each iteration, which GENCAN's
            # argument-less callback doesn't provide.
            gencan!(
                fg_closure, x;
                callback=() -> packmol_callback(nothing, cl_system, tol, iprint,
                    packmol_system.tolerance_precision, packmol_system.constraint_precision, progress_meter;
                    stall_detector=nothing, stall_tolerance=packmol_system.stall_tolerance,
                ),
                nitmax=maxit,
                nfevalmax=10 * maxit,
            )
        end

        # Update molecule positions with optimized values
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(free_mol_indices)
            packmol_system.molecule_positions[imol] = x_mol[k]
        end

        # Statistics: compute improvement of this loop relative to flast
        fx = optresult.f
        dmin = min(cl_system.fg.dmin, cl_system.cutoff)
        fimp_last = flast > zero(T) ? clamp(-100 * (fx - flast) / flast, T(-99.99), T(99.99)) : T(100)
        fimprov = bestf < typemax(T) ? clamp(-100 * (fx - bestf) / bestf, T(-99.99), T(99.99)) : T(100)
        improved = fx < bestf
        if improved
            bestf = fx
            copyto!(best_positions, packmol_system.molecule_positions)
        end
        flast = fx
        max_const = cl_system.fg.max_constraint_penalty

        finish!(progress_meter)
        @printf("\n  Function value from last loop: f = %10.5e\n", fx)
        @printf("  Best function value before: f = %10.5e\n", bestf)
        @printf("  Improvement from best function value: %8.2f %%\n", fimprov)
        @printf("  Improvement from last loop: %8.2f %%\n", fimp_last)
        @printf("  Maximum violation of target distance: %12.6f\n", tol - dmin)
        @printf("  Maximum violation of the constraints: %10.5e\n", max_const)

        # Check convergence: both tolerance and constraint precisions must be satisfied
        tol_ok = tol - dmin < packmol_system.tolerance_precision
        const_ok = max_const < packmol_system.constraint_precision
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
            println(" Current solution written to file: ", packmol_system.output_file)
        end

        # Move bad molecules once the tolerance has been fully tightened
        # (radscale == 1, matching Fortran) and this loop did not improve
        # f by at least 10%.
        if radscale == one(T) && fimp_last < T(10)
            cm_min, cm_max = compute_cm_bounds(packmol_system)
            nmoved = movebad!(
                packmol_system, cl_system.fg.fmol, free_mol_indices, mol_structure_type, RNG;
                movefrac, precision,
                cm_lo_type=cm_min, cm_hi_type=cm_max,
            )
            if nmoved > 0
                println(" Moved $nmoved bad molecules randomly to new positions.")
            end
        end
        # Re-pack optimizer variables from (possibly moved) molecule positions
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(free_mol_indices)
            x_mol[k] = packmol_system.molecule_positions[imol]
        end

        # Loosen-tolerance heuristic decay: relax radscale toward 1.0 once
        # improvement stalls (either already within tolerance precision and
        # improving less than 10%, or improving less than 2% regardless).
        # This sets the radscale movebad above will see next loop.
        if radscale > one(T)
            if (tol_ok && fimp_last < T(10)) || fimp_last < T(2)
                radscale = max(T(0.9) * radscale, one(T))
            end
        end
    end

    if !converged
        println()
        println(hash_line)
        @printf(" WARNING: packing did not converge after %d loops (best f = %.4e)\n", nloop, bestf)
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
        println(" Solution written to file: ", packmol_system.output_file)
    end

    println()
    println(dash_line)
    tend = time()
    @printf(" Running time: %12.4f seconds.\n", tend - tstart)
    println(dash_line)
    println()

    return nothing
end

#
# Aesthetic line constants (matching Fortran Packmol)
#
const dash_line = repeat('-', 80)
const hash_line = repeat('#', 80)

#
# Tracks whether `f` has stopped meaningfully improving over several
# *consecutive* SPGBox iterations within a single packing-loop chunk — a
# more robust stall signal than an instantaneous gradient/function ratio
# (tried and reverted: on this coupled, many-body objective, the projected
# gradient can look artificially small for the first iteration or two after
# a chunk starts, e.g. right after movebad! scrambles some positions, well
# before the spectral step-size estimate has calibrated — that caused chunks
# to bail out almost immediately, before making any real progress).
#
# `f_prev = nothing` marks "no baseline yet" (start of a chunk), so the
# first call is never counted as a plateau. `plateau_count` counts
# consecutive iterations whose relative improvement over the previous one
# falls below `rel_tol`; it resets to 0 on any iteration that improves by
# more than that. `is_stalled!` reports true once `plateau_count` reaches
# `max_plateau`.
#
mutable struct StallDetector{T}
    f_prev::Union{Nothing,T}
    plateau_count::Int
end
StallDetector{T}() where {T} = StallDetector{T}(nothing, 0)

function is_stalled!(detector::StallDetector{T}, f::T; rel_tol::T, max_plateau::Int=5) where {T}
    f_prev = detector.f_prev
    detector.f_prev = f
    isnothing(f_prev) && return false
    rel_improvement = f_prev > zero(T) ? (f_prev - f) / f_prev : zero(T)
    if rel_improvement < rel_tol
        detector.plateau_count += 1
    else
        detector.plateau_count = 0
    end
    return detector.plateau_count >= max_plateau
end

#
# SPGBox callback: print progress and check convergence
#
function packmol_callback(
    spgresult, cl_system, tol, iprint, tolerance_precision, constraint_precision, progress_meter;
    stall_detector::Union{Nothing,StallDetector} = nothing,
    stall_tolerance = 1e-2 * tolerance_precision,
)
    next!(progress_meter; showvalues = [
        (" Minimum distance: ", min(cl_system.fg.dmin, cl_system.cutoff)),
        (" Maximum constraint violation: ", cl_system.fg.max_constraint_penalty),
    ])
    dmin = min(cl_system.fg.dmin, cl_system.cutoff)
    tol_ok = tol - dmin < tolerance_precision
    const_ok = cl_system.fg.max_constraint_penalty < constraint_precision
    if tol_ok && const_ok
        return true
    end
    # Cut this chunk short once f has plateaued for several consecutive
    # iterations: a handful of molecules genuinely stuck (deep inside a
    # fixed structure, or wedged against others with nowhere left to go)
    # produce exactly this signature, while the bulk of the system is
    # otherwise fine — grinding through the rest of this chunk's maxit
    # budget on them cannot help; movebad! (once this chunk returns to the
    # outer loop) is what actually relocates them.
    if !isnothing(stall_detector) && is_stalled!(stall_detector, spgresult.f; rel_tol=stall_tolerance)
        return true
    end
    return false
end

@testitem "water box packing" begin
    using Packmol
    using PDBTools: read_pdb
    using StaticArrays
    using LinearAlgebra: norm

    input_file = joinpath(Packmol.src_dir, "..", "test", "input_files", "water_box_small.inp")
    sys = Packmol.read_packmol_input(input_file)
    @test sys.nmols == 100
    @test length(sys.atoms) == 300
    @test sys.tolerance == 2.0

    Packmol.packmol(sys; nloop=200, maxit=20, iprint=10)

    # Compute final atom positions
    atom_positions = Vector{SVector{3,Float64}}(undef, length(sys.atoms))
    Packmol.compute_atom_positions!(atom_positions, sys.molecule_positions, sys)

    # Check that all atoms satisfy the box constraint (penalty ≈ 0)
    st = sys.structure_types[1]
    c = st.constraints[1]  # InsideBox
    for pos in atom_positions
        @test Packmol.constraint_penalty(c, pos) ≈ 0.0 atol = 0.1
    end

    # Check minimum inter-molecular distance is close to tolerance
    natoms_per_mol = st.natoms
    dmin = let dmin = Inf
        for i in 1:length(atom_positions)
            mol_i = (i - 1) ÷ natoms_per_mol + 1
            for j in i+1:length(atom_positions)
                mol_j = (j - 1) ÷ natoms_per_mol + 1
                if mol_i != mol_j
                    d = norm(atom_positions[i] - atom_positions[j])
                    dmin = min(dmin, d)
                end
            end
        end
        dmin
    end
    @test dmin >= sys.tolerance - 0.1

    # Write output and verify
    output_file = Packmol.write_output(sys)
    @test isfile(output_file)
    output_atoms = read_pdb(output_file)
    @test length(output_atoms) == 300
    rm(output_file; force=true)
end

@testitem "water box packing (gencan)" begin
    using Packmol

    if !Packmol.gencan_gfortran_available()
        @warn "gfortran not available: skipping gencan integration test"
        @test_skip false
    else
        using PDBTools: read_pdb
        using StaticArrays
        using LinearAlgebra: norm

        input_file = joinpath(Packmol.src_dir, "..", "test", "input_files", "water_box_small.inp")
        sys = Packmol.read_packmol_input(input_file)
        @test sys.nmols == 100
        @test length(sys.atoms) == 300
        @test sys.tolerance == 2.0

        Packmol.packmol(sys; nloop=200, maxit=20, iprint=10, optimizer=:gencan)

        # Compute final atom positions
        atom_positions = Vector{SVector{3,Float64}}(undef, length(sys.atoms))
        Packmol.compute_atom_positions!(atom_positions, sys.molecule_positions, sys)

        # Check that all atoms satisfy the box constraint (penalty ≈ 0)
        st = sys.structure_types[1]
        c = st.constraints[1]  # InsideBox
        for pos in atom_positions
            @test Packmol.constraint_penalty(c, pos) ≈ 0.0 atol = 0.1
        end

        # Check minimum inter-molecular distance is close to tolerance
        natoms_per_mol = st.natoms
        dmin = let dmin = Inf
            for i in 1:length(atom_positions)
                mol_i = (i - 1) ÷ natoms_per_mol + 1
                for j in i+1:length(atom_positions)
                    mol_j = (j - 1) ÷ natoms_per_mol + 1
                    if mol_i != mol_j
                        d = norm(atom_positions[i] - atom_positions[j])
                        dmin = min(dmin, d)
                    end
                end
            end
            dmin
        end
        @test dmin >= sys.tolerance - 0.1

        # Write output and verify
        output_file = Packmol.write_output(sys)
        @test isfile(output_file)
        output_atoms = read_pdb(output_file)
        @test length(output_atoms) == 300
        rm(output_file; force=true)
    end
end