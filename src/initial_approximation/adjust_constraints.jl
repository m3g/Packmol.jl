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
    #
    # Unlike reinitialize_with_bounds.jl (which queries one trial molecule at a
    # time and can cheaply reject far-away ones before touching CellListMap at
    # all), molecule_badness.jl below does a single pairwise! call over *every*
    # free atom at once — so this system's box genuinely needs to cover the
    # whole free-molecule placement region, not just the fixed atoms. Without
    # an explicit unitcell, CellListMap treats it as non-periodic and, on the
    # very first update! (all free-atom positions at once), has to grow the
    # box from "just the fixed atoms" (its construction-time extent) out to
    # "fixed atoms ∪ every free atom" in one shot — an expensive rebuild of a
    # cell list sized however large that turns out to be. Sizing the box up
    # front from the known CM bounds avoids that surprise growth/rebuild.
    fixed_sys = if domovebad
        fixed_unitcell = if !isnothing(packmol_system.unitcell)
            packmol_system.unitcell
        else
            lo, hi = if !isnothing(cm_lo_type) && !isnothing(cm_hi_type)
                cm_lo_type, cm_hi_type
            else
                compute_cm_bounds(packmol_system)
            end
            region_lo = reduce((a, b) -> min.(a, b), lo)
            region_hi = reduce((a, b) -> max.(a, b), hi)
            fixed_atom_positions = _collect_fixed_positions(packmol_system)
            if !isempty(fixed_atom_positions)
                flo, fhi = compute_bounding_box(fixed_atom_positions)
                region_lo = min.(region_lo, flo)
                region_hi = max.(region_hi, fhi)
            end
            max_extent = zero(T)
            for st in packmol_system.structure_types
                st.fixed.fixed && continue
                for r in st.reference_coordinates
                    max_extent = max(max_extent, norm(r))
                end
            end
            margin = max_extent + packmol_system.tolerance
            T(1.2) * ((region_hi .+ margin) - (region_lo .- margin))
        end
        # Cell count scales as box volume / cutoff^D. This box's volume is set
        # by where free molecules *could* be (the whole placement region),
        # which for a large, sparse system (fixed structure much smaller than
        # the region free molecules are packed into) can imply a cell count
        # many orders of magnitude above the actual atom count — the cell
        # array itself is allocated at that size regardless of parallel vs
        # serial construction, which is what actually blows up memory here.
        # Cap it at (fixed + free) atom count, matching packmol_main.jl's own
        # cl_system cutoff cap.
        volume = fixed_unitcell isa AbstractMatrix ? abs(det(fixed_unitcell)) : prod(fixed_unitcell)
        max_particles = natoms
        cutoff = _capped_cutoff(volume, packmol_system.tolerance, max_particles, D)
        # parallel=false: CellListMap's *parallel* cell-list construction (not
        # the pairwise! computation itself) allocates per-batch auxiliary
        # arrays proportional to nthreads * ncells, which dominates memory at
        # the cell counts this box implies (millions, for a large system) —
        # confirmed 8x more memory than parallel=false for an equivalent
        # standalone system, with no speed benefit at the loop-count adjust_constraints!
        # actually runs (it typically converges in a handful of iterations).
        build_fixed_particle_system(packmol_system; unitcell=fixed_unitcell, parallel=false, cutoff)
    else
        nothing
    end
    tol = packmol_system.tolerance

    # Avoid fixed-atom overlaps during this placement phase against the same
    # (loose) distance the main packing loop's own first iteration actually
    # enforces — `radscale * tolerance` (Fortran Packmol's `discale`; see
    # `packmol_system.radscale`'s docstring) — rather than the bare
    # tolerance. Placing free molecules only `tolerance` away from the fixed
    # structure still leaves them violating that first, loosest loop, which
    # just pushes the problem into the (slower) main optimization instead of
    # avoiding it here where it's cheap.
    overlap_tol = packmol_system.radscale * tol

    # Separate, single-molecule-at-a-time overlap-check system for the
    # movebad loop below: randomize_molecule! draws one trial molecule at a
    # time (not a bulk `pairwise!` over every free atom, unlike `fixed_sys`
    # above), so it needs the `nmols=1` kind of system `overlaps_fixed`
    # expects — see `_build_overlap_check_system`. Without this, a molecule
    # randomized away from one bad position had no way to avoid landing on
    # top of the fixed structure again, which for a large/dense fixed
    # structure (e.g. a real protein) could take many loops to resolve by
    # chance, or never resolve within `nloop` at all.
    randomize_fixed_sys, randomize_fixed_lo, randomize_fixed_hi = domovebad ?
        _build_overlap_check_system(packmol_system, overlap_tol) :
        (nothing, zero(SVector{D,T}), zero(SVector{D,T}))

    # Build molecule → structure type / first-atom-index mappings
    mol_structure_type = _build_mol_structure_type(packmol_system)
    mol_iat_first = _build_mol_iat_first(packmol_system)

    # Persistent, cumulative per-molecule badness cache. Molecules are
    # geometrically independent in this constraint-only phase (they interact
    # only through their own constraints and overlap with *fixed* atoms, never
    # with each other), so once a molecule's badness is known it stays valid
    # until that molecule itself is moved — no need to recompute it on every
    # loop just because *other* molecules are still being adjusted.
    badness = isnothing(buffers) ? zeros(T, packmol_system.nmols) : buffers.badness
    fill!(badness, zero(T))

    # Scratch buffers for the fixed-atom overlap check, sized to the largest
    # possible active set (all free molecules, on loop 1) and reused (as a
    # prefix) on every subsequent, smaller loop.
    n_free_atoms_max = sum(
        st.fixed.fixed ? 0 : st.number_of_molecules * st.natoms
        for st in packmol_system.structure_types
    )
    active_free_atoms = Vector{SVector{D,T}}(undef, n_free_atoms_max)
    active_free_atom_mol = Vector{Int}(undef, n_free_atoms_max)

    # `fixed_sys` (built above only when domovebad=true, and only if there
    # are fixed atoms and avoid_overlap is on) is what makes badness/"bad
    # molecules" below include fixed-structure overlap on top of geometric
    # constraint violations — not just a cosmetic label: it's why this same
    # function's two calls in set_initial_approximation! (domovebad=false,
    # then domovebad=true) can report very different "bad" counts for what
    # looks like the same problem. Say so up front so the two phases don't
    # read as inconsistent.
    checking_str = isnothing(fixed_sys) ? "constraints" : "constraints + fixed overlap"
    println("  Adjusting initial point to fit $checking_str")

    # The active set starts as every free molecule (all need their first
    # evaluation); it only shrinks from here — a molecule leaves it as soon
    # as it's satisfied, and only re-enters if explicitly randomized.
    active_mols = copy(free_mol_indices)

    for iloop in 1:nloop
        nactive = length(active_mols)

        # Set up optimizer variables from current molecule positions. Reuse
        # the pre-allocated buffer only when the active set is still at full
        # size (loop 1, typically the one expensive pass); once it has
        # shrunk, allocate a right-sized one instead — cheap, since by then
        # nactive is small.
        x = (!isnothing(buffers) && nactive == nfree) ? buffers.x : Vector{T}(undef, nactive * 2 * D)
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(active_mols)
            x_mol[k] = packmol_system.molecule_positions[imol]
        end
        auxvecs = (!isnothing(buffers) && nactive == nfree) ? buffers.vaux : SPGBox.VAux(x, zero(T))

        # Run a short constraint-only optimization, touching only active_mols
        spgresult = spgbox!(
            (g, x) -> constraint_only_fg_for_mols!(
                g, x, fg_output, packmol_system, atom_positions, active_mols, mol_structure_type, mol_iat_first,
            ),
            x;
            vaux=auxvecs,
            nitmax=opt_nit,
            nfevalmax=10 * opt_nit,
            callback=(result) -> result.f < precision,
        )

        # Update molecule positions from optimizer
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(active_mols)
            packmol_system.molecule_positions[imol] = x_mol[k]
        end

        # Compute atom positions for per-molecule evaluation, only for active_mols
        compute_positions_and_constraints_for_mols!(
            atom_positions, fg_output, packmol_system, active_mols, mol_structure_type, mol_iat_first,
        )

        # Update the badness cache for active_mols only; entries for molecules
        # outside active_mols are untouched — they haven't moved, so their
        # last-computed badness is still correct.
        molecule_badness_for_mols!(
            badness, packmol_system, atom_positions, fixed_sys, overlap_tol,
            active_mols, mol_structure_type, mol_iat_first, active_free_atoms, active_free_atom_mol,
        )

        # Check if all constraints are satisfied and no overlaps with fixed
        total_badness = sum(badness[imol] for imol in free_mol_indices)
        if total_badness < precision
            @printf("  Adjustment for %s converged at loop %d (f = %.2e)\n", checking_str, iloop, total_badness)
            return packmol_system
        end

        # Identify bad molecules — only among active_mols: anything outside
        # it is known-good and untouched this loop.
        bad_mols = Int[]
        bad_vals = T[]
        for imol in active_mols
            if badness[imol] > precision
                push!(bad_mols, imol)
                push!(bad_vals, badness[imol])
            end
        end
        nbad = length(bad_mols)

        if nbad == 0
            @printf("  All molecules satisfy %s at loop %d\n", checking_str, iloop)
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
                    fixed_sys = randomize_fixed_sys,
                    fixed_lo = randomize_fixed_lo,
                    fixed_hi = randomize_fixed_hi,
                    tol = overlap_tol,
                )
            else
                randomize_molecule!(packmol_system, bad_mols[i], st, RNG;
                    fixed_sys = randomize_fixed_sys,
                    fixed_lo = randomize_fixed_lo,
                    fixed_hi = randomize_fixed_hi,
                    tol = overlap_tol,
                )
            end
        end

        # Next loop's active set: every molecule still in violation (whether
        # just randomized or left as-is because movefrac capped how many
        # could move) — everything else is now known-good and stays excluded.
        active_mols = bad_mols
    end

    # Final evaluation (full system, for an accurate warning message). fixed_sys's
    # output may have been shrunk by molecule_badness_for_mols!'s active-set
    # resizing above; molecule_badness (unrestricted) indexes it by global
    # molecule id, so it must be restored to its full nmols size first.
    if !isnothing(fixed_sys)
        CellListMap.resize_output!(fixed_sys, packmol_system.nmols)
    end
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)
    badness = molecule_badness(packmol_system, atom_positions, fixed_sys, tol; buffers)
    total_badness = sum(badness[imol] for imol in free_mol_indices)
    @printf("  WARNING: adjustment for %s did not fully converge after %d loops (f = %.2e)\n", checking_str, nloop, total_badness)
    return packmol_system
end
