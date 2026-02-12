#
# Functions related to the generation of the initial approximation for the packing problem.
#

#
# Bounding box computed from actual atom positions.
# Used after molecule placement to determine the region for CellListMap.
#
function compute_bounding_box(atom_positions::Vector{SVector{D,T}}) where {D,T}
    lo = SVector{D,T}(ntuple(_ -> typemax(T), D))
    hi = SVector{D,T}(ntuple(_ -> typemin(T), D))
    for x in atom_positions
        lo = min.(lo, x)
        hi = max.(hi, x)
    end
    return lo, hi
end

# Default sidemax: half-width of the large box for initial placement.
const DEFAULT_SIDEMAX = 1000.0

#
# Check if a molecule placed at position `mp` overlaps with any fixed atom.
# Returns true if any atom of the molecule is within `tol` of any fixed atom.
# Uses CellListMap cross-computation: fixed_sys is a ParticleSystem1 built from
# fixed atom positions. `mol_positions` is a preallocated buffer.
#
function overlaps_fixed(
    mp::MoleculePosition{D,T},
    ref_coords::Vector{SVector{D,T}},
    fixed_sys::CellListMap.ParticleSystem1,
    mol_positions::Vector{SVector{D,T}};
    itask
) where {D,T}
    R = eulermat(mp.angles)
    for (j, r) in enumerate(ref_coords)
        mol_positions[j] = R * r + mp.cm
    end
    n = length(ref_coords)
    fixed_sys.output[itask] = zero(eltype(fixed_sys.output))
    overlap = pairwise!(
        (pair, overlap) -> begin
            overlap[itask] += one(eltype(overlap))
            return overlap
        end,
        @view(mol_positions[1:n]), fixed_sys;
        update_lists=false,      # fixed cell lists don't change
    )
    return overlap[itask] > 0
end

# Fallback when there are no fixed atoms (no ParticleSystem built)
function overlaps_fixed(
    mp::MoleculePosition{D,T},
    ref_coords::Vector{SVector{D,T}},
    ::Nothing,
    mol_positions::Vector{SVector{D,T}};
    itask
) where {D,T}
    return false
end

#
# Check if a molecule placement satisfies all constraints of its structure type.
# Returns true if all atom positions have zero constraint penalty.
#
function satisfies_constraints(
    mp::MoleculePosition{D,T},
    st::StructureType{D,T},
    mol_positions::Vector{SVector{D,T}},
) where {D,T}
    isempty(st.constraints) && return true
    R = eulermat(mp.angles)
    for (j, r) in enumerate(st.reference_coordinates)
        x = R * r + mp.cm
        mol_positions[j] = x
        for ic in st.atom_constraints[j]
            c = st.constraints[ic]
            if constraint_penalty(c, x) > zero(T)
                return false
            end
        end
    end
    return true
end

#
# Total constraint penalty for a molecule placement.
# Returns the sum of constraint penalties for all atoms of the molecule.
#
function constraint_penalty_sum(
    mp::MoleculePosition{D,T},
    st::StructureType{D,T},
) where {D,T}
    isempty(st.constraints) && return zero(T)
    R = eulermat(mp.angles)
    fpen = zero(T)
    for (j, r) in enumerate(st.reference_coordinates)
        x = R * r + mp.cm
        for ic in st.atom_constraints[j]
            c = st.constraints[ic]
            fpen += constraint_penalty(c, x)
        end
    end
    return fpen
end

#
# Compute per-structure-type center-of-mass bounds from current molecule positions.
# Returns (cm_min, cm_max) vectors indexed by structure type.
# Used after constraint fitting to determine where molecules should be placed.
#
function compute_cm_bounds(packmol_system::PackmolSystem{D,T}) where {D,T}
    ntypes = length(packmol_system.structure_types)
    cm_min = [SVector{D,T}(ntuple(_ -> typemax(T), D)) for _ in 1:ntypes]
    cm_max = [SVector{D,T}(ntuple(_ -> typemin(T), D)) for _ in 1:ntypes]
    imol = 0
    for (ist, st) in enumerate(packmol_system.structure_types)
        for _ in 1:st.number_of_molecules
            imol += 1
            if !st.fixed.fixed
                cm = packmol_system.molecule_positions[imol].cm
                cm_min[ist] = min.(cm_min[ist], cm)
                cm_max[ist] = max.(cm_max[ist], cm)
            end
        end
    end
    return cm_min, cm_max
end

#
# Estimate per-structure-type CM bounds by evaluating constraint penalty + overlap
# with fixed atoms at the current (random) molecule positions and keeping
# the best fraction of molecules per type.
#
# Returns (cm_min, cm_max) vectors indexed by structure type, computed from
# the CMs of the `best_fraction` molecules with lowest badness per type.
#
function estimate_cm_bounds(
    packmol_system::PackmolSystem{D,T};
    best_fraction::T=T(0.1),
) where {D,T}
    # Compute atom positions and per-molecule badness
    natoms = length(packmol_system.atoms)
    atom_positions = Vector{SVector{D,T}}(undef, natoms)
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)
    fixed_sys = build_fixed_particle_system(packmol_system)
    tol = packmol_system.tolerance
    badness = molecule_badness(packmol_system, atom_positions, fixed_sys, tol)

    ntypes = length(packmol_system.structure_types)
    cm_min = [SVector{D,T}(ntuple(_ -> typemax(T), D)) for _ in 1:ntypes]
    cm_max = [SVector{D,T}(ntuple(_ -> typemin(T), D)) for _ in 1:ntypes]

    imol = 0
    for (ist, st) in enumerate(packmol_system.structure_types)
        if st.fixed.fixed
            imol += st.number_of_molecules
            continue
        end
        # Collect (molecule_index, badness) for this type
        mol_indices = Vector{Int}(undef, st.number_of_molecules)
        mol_badness = Vector{T}(undef, st.number_of_molecules)
        for i in 1:st.number_of_molecules
            imol += 1
            mol_indices[i] = imol
            mol_badness[i] = badness[imol]
        end
        # Sort by badness (ascending) and take best fraction
        perm = sortperm(mol_badness)
        nbest = max(1, round(Int, best_fraction * st.number_of_molecules))
        for k in 1:nbest
            cm = packmol_system.molecule_positions[mol_indices[perm[k]]].cm
            cm_min[ist] = min.(cm_min[ist], cm)
            cm_max[ist] = max.(cm_max[ist], cm)
        end
    end
    return cm_min, cm_max
end

#
# Initialize molecule positions randomly within constraint bounding box
# (or within the unit cell for PBC), and center reference coordinates
# at the origin (required for chain rule).
#
# This is just the first rough placement — no constraint checking is done here.
# The subsequent steps will estimate bounds from the best molecules and
# re-initialize within those bounds.
#
function initialize_molecules!(packmol_system::PackmolSystem{D,T}, RNG) where {D,T}
    # Center reference coordinates at origin (required for chain rule)
    for st in packmol_system.structure_types
        if !st.fixed.fixed
            cm = mean(st.reference_coordinates)
            st.reference_coordinates .-= Ref(cm)
        end
    end

    # Determine the placement region
    has_pbc = !isnothing(packmol_system.unitcell)
    if has_pbc
        uc = packmol_system.unitcell
        center = packmol_system.unitcell_center
    end

    # Initial placement: use a large box centered at origin (following Fortran
    # Packmol's sidemax approach). The constraint-only optimization that follows
    # will move molecules to feasible positions regardless of constraint type.
    sidemax = T(DEFAULT_SIDEMAX)

    # Compute molecule index offset for each structure type so threads
    # can determine the correct slot without a shared counter.
    imol_offset = 0
    @sync for st in packmol_system.structure_types
        st_offset = imol_offset
        for (_, irange) in enumerate(chunks(1:st.number_of_molecules; n=Threads.nthreads()))
            task_seed = rand(RNG, UInt64)
            @spawn begin
                task_rng = typeof(RNG)(task_seed)
                for i in irange
                    imol_local = st_offset + i
                    if st.fixed.fixed
                        packmol_system.molecule_positions[imol_local] = st.fixed.position
                    else
                        if has_pbc
                            frac = SVector{D,T}(ntuple(_ -> rand(task_rng, T) - T(0.5), D))
                            cm = SVector{D,T}(uc * frac) + center
                        else
                            cm = SVector{D,T}(ntuple(_ -> sidemax * (T(2) * rand(task_rng, T) - one(T)), D))
                        end
                        angles = SVector{D,T}(ntuple(_ -> T(2π) * rand(task_rng, T), D))
                        packmol_system.molecule_positions[imol_local] = MoleculePosition(cm, angles)
                    end
                end
            end
        end
        imol_offset += st.number_of_molecules
    end
    return packmol_system
end

#
# Collect Cartesian positions of all fixed atoms.
# Returns an empty vector if avoid_overlap is false or there are no fixed atoms.
#
function _collect_fixed_positions(packmol_system::PackmolSystem{D,T}) where {D,T}
    if !packmol_system.avoid_overlap
        return SVector{D,T}[]
    end
    fixed_atom_positions = SVector{D,T}[]
    for st in packmol_system.structure_types
        if st.fixed.fixed
            mp = st.fixed.position
            R = eulermat(mp.angles)
            for r in st.reference_coordinates
                push!(fixed_atom_positions, R * r + mp.cm)
            end
        end
    end
    return fixed_atom_positions
end

#
# Build a mapping from molecule index to structure type index.
#
function _build_mol_structure_type(packmol_system::PackmolSystem)
    mol_structure_type = Vector{Int}(undef, packmol_system.nmols)
    imol = 0
    for (ist, st) in enumerate(packmol_system.structure_types)
        for _ in 1:st.number_of_molecules
            imol += 1
            mol_structure_type[imol] = ist
        end
    end
    return mol_structure_type
end

#
# Build a ParticleSystem1 from pre-collected fixed atom positions.
# Returns `nothing` if positions are empty.
#
function _build_fixed_particle_system(fixed_atom_positions::Vector{SVector{D,T}}, tol::T, ::Type{T}; nthreads=Threads.nthreads()) where {D,T}
    isempty(fixed_atom_positions) && return nothing
    return ParticleSystem(
        xpositions=fixed_atom_positions,
        cutoff=tol,
        output=zeros(T, nthreads),
        parallel=false,
    )
end

#
# Build a ParticleSystem1 from fixed atom positions for efficient overlap checks.
# Returns `nothing` if there are no fixed atoms or avoid_overlap is false.
#
function build_fixed_particle_system(packmol_system::PackmolSystem{D,T}) where {D,T}
    fixed_atom_positions = _collect_fixed_positions(packmol_system)
    return _build_fixed_particle_system(fixed_atom_positions, packmol_system.tolerance, T)
end

#
# Compute per-molecule "badness" = constraint penalty + overlap with fixed molecules.
# Only considers constraint violations and distance violations against fixed atoms.
# Uses CellListMap cross-computation for efficient overlap detection.
# Returns a vector of badness values indexed by molecule index.
#
function molecule_badness(
    packmol_system::PackmolSystem{D,T},
    atom_positions::Vector{SVector{D,T}},
    fixed_sys::Union{Nothing,CellListMap.ParticleSystem1},
    tol::T,
) where {D,T}
    has_pbc = !isnothing(packmol_system.unitcell)
    nmols = packmol_system.nmols
    badness = zeros(T, nmols)

    ntasks = isnothing(fixed_sys) ? Threads.nthreads() : length(fixed_sys.output)
    max_natoms = maximum(st.natoms for st in packmol_system.structure_types; init=0)

    # Compute constraint penalties and fixed-atom overlaps per molecule
    imol_offset = 0
    atom_offset = 0
    @sync for st in packmol_system.structure_types
        st_mol_offset = imol_offset
        st_atom_offset = atom_offset

        for (itask, irange) in enumerate(chunks(1:st.number_of_molecules; n=ntasks))
            @spawn begin
                mol_positions = Vector{SVector{D,T}}(undef, max_natoms)
                for i in irange
                    imol = st_mol_offset + i
                    iat_first = st_atom_offset + (i - 1) * st.natoms + 1
                    for j in 1:st.natoms
                        iat = iat_first + j - 1
                        x = atom_positions[iat]
                        if has_pbc
                            x = wrap_to_center(x, packmol_system.unitcell, packmol_system.unitcell_center)
                        end
                        for ic in packmol_system.atoms[iat].constraints
                            c = st.constraints[ic]
                            badness[imol] += constraint_penalty(c, x)
                        end
                    end
                    # Per-molecule overlap with fixed atoms via CellListMap cross-computation
                    if !isnothing(fixed_sys) && !st.fixed.fixed
                        for j in 1:st.natoms
                            mol_positions[j] = atom_positions[iat_first+j-1]
                        end
                        fixed_sys.output[itask] = zero(T)
                        overlap = pairwise!(
                            (pair, output) -> begin
                                output[itask] += (pair.d - tol)^2
                                return output
                            end,
                            @view(mol_positions[1:st.natoms]), fixed_sys;
                            update_lists=false,
                        )
                        badness[imol] += overlap[itask]
                    end
                end
            end
        end
        imol_offset += st.number_of_molecules
        atom_offset += st.number_of_molecules * st.natoms
    end

    return badness
end

#
# Move the worst molecules to new random positions.
# Scans fmol for free molecules with fmol > precision, computes the number
# of bad molecules and moves a fraction of them randomly.
# No extra arrays are allocated: fmol itself carries all the information.
# Following the Fortran Packmol heuristic (heuristics.f90 movebad subroutine).
#
function movebad!(
    packmol_system::PackmolSystem{D,T},
    fmol::Vector{T},
    free_mol_indices::Vector{Int},
    mol_structure_type::Vector{Int},
    RNG;
    movefrac::T=T(0.05),
    precision::T=T(1e-2),
    cm_lo_type::Union{Nothing,Vector{SVector{D,T}}} = nothing,
    cm_hi_type::Union{Nothing,Vector{SVector{D,T}}} = nothing,
) where {D,T}
    nfree = length(free_mol_indices)
    # Count bad molecules and find fmol range among them
    nbad = 0
    fmol_max = zero(T)
    for imol in free_mol_indices
        if fmol[imol] > precision / packmol_system.nmols
            nbad += 1
            fmol_max = max(fmol_max, fmol[imol])
        end
    end
    nbad == 0 && return 0
    # Number of molecules to move
    frac = min(movefrac, nbad / nfree)
    nmove = max(1, min(nbad, round(Int, frac * nfree)))
    # Move molecules randomly: probability of moving is proportional
    # to fmol value (worse molecules are more likely to be moved).
    nmoved = 0
    for imol in free_mol_indices
        nmoved >= nmove && break
        if fmol[imol] > precision / packmol_system.nmols
            # Probability increases with fmol value: always move the worst,
            # linearly decreasing probability for better molecules
            prob = fmol[imol] / fmol_max
            if rand(RNG, T) < prob
                ist = mol_structure_type[imol]
                st = packmol_system.structure_types[ist]
                if !isnothing(cm_lo_type) && !isnothing(cm_hi_type)
                    lo = cm_lo_type[ist]
                    hi = cm_hi_type[ist]
                    has_valid_bounds = all(lo .< hi)
                    randomize_molecule!(packmol_system, imol, st, RNG;
                        cm_lo = has_valid_bounds ? lo : nothing,
                        cm_hi = has_valid_bounds ? hi : nothing,
                    )
                else
                    randomize_molecule!(packmol_system, imol, st, RNG)
                end
                nmoved += 1
            end
        end
    end
    return nmoved
end

#
# Randomly re-place a molecule within the placement region.
# Uses per-structure-type bounding box for non-PBC placement.
#
function randomize_molecule!(
    packmol_system::PackmolSystem{D,T},
    imol::Int,
    st::StructureType{D,T},
    RNG;
    cm_lo::Union{Nothing,SVector{D,T}} = nothing,
    cm_hi::Union{Nothing,SVector{D,T}} = nothing,
) where {D,T}
    has_pbc = !isnothing(packmol_system.unitcell)
    if has_pbc
        uc = packmol_system.unitcell
        center = packmol_system.unitcell_center
        frac = SVector{D,T}(ntuple(_ -> rand(RNG, T) - T(0.5), D))
        cm = SVector{D,T}(uc * frac) + center
    elseif !isnothing(cm_lo) && !isnothing(cm_hi)
        extent = cm_hi - cm_lo
        cm = cm_lo + SVector{D,T}(ntuple(d -> rand(RNG, T) * extent[d], D))
    else
        sidemax = T(DEFAULT_SIDEMAX)
        cm = SVector{D,T}(ntuple(_ -> sidemax * (T(2) * rand(RNG, T) - one(T)), D))
    end
    angles = SVector{D,T}(ntuple(_ -> T(2π) * rand(RNG, T), D))
    packmol_system.molecule_positions[imol] = MoleculePosition(cm, angles)
    return nothing
end

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

        for (itask, irange) in enumerate(chunks(1:st.number_of_molecules; n=Threads.nthreads()))
            task_seed = rand(RNG, UInt64)
            @spawn begin
                task_rng = typeof(RNG)(task_seed)
                mol_positions = Vector{SVector{D,T}}(undef, max_natoms)
                for i in irange
                    imol = st_offset + i
                    # Start fresh: the goal is to replace the current position
                    # (which may be edge-clustered from constraint optimization)
                    # with a uniformly distributed random placement.
                    best_mp = packmol_system.molecule_positions[imol]
                    best_fx = typemax(T)

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
                           overlaps_fixed(mp, st.reference_coordinates, fixed_sys, mol_positions; itask)
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

#
# Constraint-only objective function + gradient for pre-optimization.
# Evaluates only constraint penalties (no atom-atom distance penalties),
# so molecules can move freely to satisfy their geometric constraints
# before the main packing optimization.
#
function constraint_only_fg!(
    g, x,
    fg_output::InteratomicDistanceFG{D,T},
    packmol_system::PackmolSystem{D,T},
    atom_positions::Vector{SVector{D,T}},
    free_mol_indices::Vector{Int},
) where {D,T}
    # Unpack optimizer variables into free molecule slots
    x_mol = reinterpret(MoleculePosition{D,T}, x)
    for (k, imol) in enumerate(free_mol_indices)
        packmol_system.molecule_positions[imol] = x_mol[k]
    end
    # Reset fg output (no CellListMap to do it for us)
    CellListMap.reset_output!(fg_output)
    # Single pass: compute atom positions AND constraint penalties/gradients
    compute_positions_and_constraints!(atom_positions, fg_output, packmol_system.molecule_positions, packmol_system)
    # Chain rule: Cartesian → molecule DOF gradients (for all molecules)
    chain_rule!(fg_output, packmol_system)
    # Pack only free molecule gradients into optimizer gradient vector
    g_mol = reinterpret(MoleculePosition{D,T}, g)
    for (k, imol) in enumerate(free_mol_indices)
        g_mol[k] = fg_output.g[imol]
    end
    return fg_output.f
end

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
    atom_positions = Vector{SVector{D,T}}(undef, natoms)
    fg_output = InteratomicDistanceFG{D,T}(packmol_system)

    # Build CellListMap system for fixed atom positions (for overlap checks)
    fixed_sys = build_fixed_particle_system(packmol_system)
    tol = packmol_system.tolerance

    # Build molecule → structure type mapping
    mol_structure_type = _build_mol_structure_type(packmol_system)

    println("  Adjusting initial point to fit the constraints")

    for iloop in 1:nloop
        # Set up optimizer variables from current molecule positions
        x = Vector{T}(undef, nfree * 2 * D)
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(free_mol_indices)
            x_mol[k] = packmol_system.molecule_positions[imol]
        end
        auxvecs = SPGBox.VAux(x, zero(T))

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
        badness = molecule_badness(packmol_system, atom_positions, fixed_sys, tol)

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
    badness = molecule_badness(packmol_system, atom_positions, fixed_sys, tol)
    total_badness = sum(badness[imol] for imol in free_mol_indices)
    @printf("  WARNING: constraint adjustment did not fully converge after %d loops (f = %.2e)\n", nloop, total_badness)
    return packmol_system
end
