#
# Functions related to the generation of the initial approximation for the packing problem.
#

#
# Bounding box computation from constraints
#
constraint_bounds(c::Box) = (c.center - c.sides / 2, c.center + c.sides / 2)
constraint_bounds(c::Cube) = (c.center .- c.side / 2, c.center .+ c.side / 2)
constraint_bounds(c::Sphere) = (c.center .- c.radius, c.center .+ c.radius)

function compute_bounding_box(packmol_system::PackmolSystem{D,T}) where {D,T}
    lo = SVector{D,T}(ntuple(_ -> typemax(T), D))
    hi = SVector{D,T}(ntuple(_ -> typemin(T), D))
    for st in packmol_system.structure_types
        for c in st.constraints
            clo, chi = constraint_bounds(c)
            lo = min.(lo, clo)
            hi = max.(hi, chi)
        end
    end
    return lo, hi
end

#
# Check if a molecule placed at position `mp` overlaps with any fixed atom.
# Returns true if any atom of the molecule is within `tol` of any fixed atom.
#
function overlaps_fixed(
    mp::MoleculePosition{D,T},
    ref_coords::Vector{SVector{D,T}},
    fixed_atom_positions::Vector{SVector{D,T}},
    tol::T,
) where {D,T}
    isempty(fixed_atom_positions) && return false
    tol2 = tol^2
    R = eulermat(mp.angles...)
    for r in ref_coords
        pos = R * r + mp.cm
        for fp in fixed_atom_positions
            d2 = sum((pos - fp) .^ 2)
            if d2 < tol2
                return true
            end
        end
    end
    return false
end

#
# Initialize molecule positions randomly within constraint bounding box
# (or within the unit cell for PBC), and center reference coordinates
# at the origin (required for chain rule).
#
function initialize_molecules!(packmol_system::PackmolSystem{D,T}, RNG) where {D,T}
    # Center reference coordinates at origin (required for chain rule)
    for st in packmol_system.structure_types
        if !st.fixed.fixed
            cm = mean(st.reference_coordinates)
            st.reference_coordinates .-= Ref(cm)
        end
    end

    # Precompute atom positions of all fixed molecules
    fixed_atom_positions = SVector{D,T}[]
    if packmol_system.avoid_overlap
        for st in packmol_system.structure_types
            if st.fixed.fixed
                mp = st.fixed.position
                R = eulermat(mp.angles...)
                for r in st.reference_coordinates
                    push!(fixed_atom_positions, R * r + mp.cm)
                end
            end
        end
    end
    tol = packmol_system.tolerance

    # Determine the placement region
    has_pbc = !isnothing(packmol_system.unitcell)
    if has_pbc
        # For PBC, place randomly within the unit cell centered at unitcell_center
        uc = packmol_system.unitcell
        center = packmol_system.unitcell_center
    else
        # For non-PBC, place within the constraint bounding box
        lo, hi = compute_bounding_box(packmol_system)
    end
    max_attempts = 1000
    imol = 0
    for st in packmol_system.structure_types
        for _ in 1:st.number_of_molecules
            imol += 1
            if st.fixed.fixed
                packmol_system.molecule_positions[imol] = st.fixed.position
            else
                for attempt in 1:max_attempts
                    if has_pbc
                        # Random fractional coords in [-0.5, 0.5) → Cartesian, centered at unitcell_center
                        frac = SVector{D,T}(ntuple(_ -> rand(RNG, T) - T(0.5), D))
                        cm = SVector{D,T}(uc * frac) + center
                    else
                        extent = hi - lo
                        cm = lo + SVector{D,T}(ntuple(d -> rand(RNG, T) * extent[d], D))
                    end
                    angles = SVector{D,T}(ntuple(_ -> T(2π) * rand(RNG, T), D))
                    mp = MoleculePosition(cm, angles)
                    if !packmol_system.avoid_overlap || !overlaps_fixed(mp, st.reference_coordinates, fixed_atom_positions, tol)
                        packmol_system.molecule_positions[imol] = mp
                        break
                    end
                    if attempt == max_attempts
                        @warn "Could not place molecule $imol without overlap after $max_attempts attempts"
                        packmol_system.molecule_positions[imol] = mp
                    end
                end
            end
        end
    end
    return packmol_system
end

#
# Compute per-molecule "badness" = constraint penalty + overlap with fixed molecules.
# Only considers constraint violations and distance violations against fixed atoms.
# Returns a vector of badness values indexed by molecule index.
#
function molecule_badness(
    packmol_system::PackmolSystem{D,T},
    atom_positions::Vector{SVector{D,T}},
    fixed_atom_positions::Vector{SVector{D,T}},
    tol::T,
) where {D,T}
    has_pbc = !isnothing(packmol_system.unitcell)
    nmols = packmol_system.nmols
    badness = zeros(T, nmols)
    tol2 = tol^2

    # Compute constraint penalties and fixed-molecule overlaps per molecule
    iat = 0
    imol = 0
    for st in packmol_system.structure_types
        for _ in 1:st.number_of_molecules
            imol += 1
            for j in 1:st.natoms
                iat += 1
                x = atom_positions[iat]
                if has_pbc
                    x = wrap_to_center(x, packmol_system.unitcell, packmol_system.unitcell_center)
                end
                # Constraint penalty
                for ic in packmol_system.atoms[iat].constraints
                    c = st.constraints[ic]
                    badness[imol] += constraint_penalty(c, x)
                end
                # Overlap with fixed atoms
                for fp in fixed_atom_positions
                    d2 = sum((atom_positions[iat] - fp) .^ 2)
                    if d2 < tol2
                        badness[imol] += (sqrt(d2) - tol)^2
                    end
                end
            end
        end
    end

    return badness
end

#
# Randomly re-place a molecule within the placement region.
#
function randomize_molecule!(
    packmol_system::PackmolSystem{D,T},
    imol::Int,
    RNG,
    lo::SVector{D,T},
    hi::SVector{D,T},
    has_pbc::Bool,
) where {D,T}
    if has_pbc
        uc = packmol_system.unitcell
        center = packmol_system.unitcell_center
        frac = SVector{D,T}(ntuple(_ -> rand(RNG, T) - T(0.5), D))
        cm = SVector{D,T}(uc * frac) + center
    else
        extent = hi - lo
        cm = lo + SVector{D,T}(ntuple(d -> rand(RNG, T) * extent[d], D))
    end
    angles = SVector{D,T}(ntuple(_ -> T(2π) * rand(RNG, T), D))
    packmol_system.molecule_positions[imol] = MoleculePosition(cm, angles)
    return nothing
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
    # Compute Cartesian atomic coordinates from ALL molecule DOFs
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)
    # Reset fg output (no CellListMap to do it for us)
    CellListMap.reset_output!(fg_output)
    # Add constraint penalties and gradients
    constraint_fg!(fg_output, atom_positions, packmol_system)
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

    # Precompute fixed atom positions for overlap checks
    fixed_atom_positions = SVector{D,T}[]
    if packmol_system.avoid_overlap
        for st in packmol_system.structure_types
            if st.fixed.fixed
                mp = st.fixed.position
                R = eulermat(mp.angles...)
                for r in st.reference_coordinates
                    push!(fixed_atom_positions, R * r + mp.cm)
                end
            end
        end
    end
    tol = packmol_system.tolerance

    # Determine placement region for randomization
    has_pbc = !isnothing(packmol_system.unitcell)
    if has_pbc
        lo = hi = zero(SVector{D,T}) # unused with PBC
    else
        lo, hi = compute_bounding_box(packmol_system)
    end

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
        badness = molecule_badness(packmol_system, atom_positions, fixed_atom_positions, tol)

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
        nmove = max(1, min(nbad, round(Int, movefrac * nfree)))
        if iloop % iprint == 0 || iloop == 1
            @printf("  Loop %4d: f = %10.4e  bad molecules: %d/%d  moving: %d\n",
                iloop, spgresult.f, nbad, nfree, nmove)
        end
        for i in 1:nmove
            randomize_molecule!(packmol_system, bad_mols[i], RNG, lo, hi, has_pbc)
        end
    end

    # Final evaluation
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)
    badness = molecule_badness(packmol_system, atom_positions, fixed_atom_positions, tol)
    total_badness = sum(badness[imol] for imol in free_mol_indices)
    @printf("  WARNING: constraint adjustment did not fully converge after %d loops (f = %.2e)\n", nloop, total_badness)
    return packmol_system
end
