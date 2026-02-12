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
    maxit::Int=100,
    movefrac::T=T(0.05),
    seed::Int=packmol_system.seed,
    restart::Bool=false,
) where {D,T}
    # Initialize RNG and molecule positions
    RNG = Random.Xoshiro(seed)

    tstart = time()    

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

    # Pre-optimization: move molecules to satisfy geometric constraints
    # before the main packing (no distance penalties)
    if restart
        initialize_molecules!(packmol_system, RNG)
    else
        # Step 1: Random positions in large sidemax box (done by initialize_molecules!)
        initialize_molecules!(packmol_system, RNG)
        # Step 2: Constraint-only optimization (no movebad) to push molecules
        # toward feasible regions, then compute CM bounds from the result
        adjust_constraints!(packmol_system, free_mol_indices, RNG; domovebad=false)
        cm_min, cm_max = compute_cm_bounds(packmol_system)
        # Step 3: Re-initialize molecules within the estimated bounds
        # (best-of-N tries per molecule, checking constraints + fixed overlap)
        reinitialize_with_bounds!(packmol_system, free_mol_indices, RNG;
            precision=packmol_system.constraint_precision,
            cm_min=cm_min, cm_max=cm_max,
        )
        # Step 4: Full constraint adjustment with the solver and movebad,
        # randomizing bad molecules within the CM bounds from step 2
        if packmol_system.adjust_constraints_on_init
            adjust_constraints!(packmol_system, free_mol_indices, RNG;
                cm_lo_type=cm_min, cm_hi_type=cm_max,
            )
        end
    end

    # check mode: write the initial approximation and return
    if packmol_system.check
        println("Check mode: writing initial approximation to $(packmol_system.output_file)")
        if !isempty(packmol_system.output_file)
            write_output(packmol_system)
        end
        return nothing
    end

    # Compute initial atom positions
    natoms = length(packmol_system.atoms)
    atom_positions = Vector{SVector{D,T}}(undef, natoms)
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)

    # Determine unit cell for CellListMap
    has_pbc = !isnothing(packmol_system.unitcell)
    if has_pbc
        # PBC mode: use the actual unit cell
        unitcell = packmol_system.unitcell
    else
        # Non-PBC mode: inflate bounding box so CellListMap treats it as a large box
        lo, hi = compute_bounding_box(atom_positions)
        box_size = hi - lo
        unitcell = T(3) * box_size
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
    ncells = min(volume / packing_tol^D, natoms)
    #cutoff = max(packing_tol, (volume / ncells)^(one(T) / D))
    cutoff = packing_tol

    # Set up CellListMap
    fg_output = InteratomicDistanceFG{D,T}(packmol_system)
    println("Setting up cell lists ($natoms atoms, cutoff = $(@sprintf("%.2f", cutoff)))...")
    cl_system = ParticleSystem(
        xpositions=atom_positions,
        unitcell=unitcell,
        cutoff=cutoff,
        output=fg_output,
        output_name=:fg,
        parallel=parallel,
    )

    # Set up optimization variables: only free molecule DOFs
    x = Vector{T}(undef, nfree * 2 * D)
    x_mol = reinterpret(MoleculePosition{D,T}, x)
    for (k, imol) in enumerate(free_mol_indices)
        x_mol[k] = packmol_system.molecule_positions[imol]
    end
    auxvecs = SPGBox.VAux(x, zero(T))

    # Placement region for movebad! randomization
    mol_structure_type = _build_mol_structure_type(packmol_system)
    precision = packmol_system.tolerance_precision

    # Outer packing loop (following Fortran Packmol gencanloop):
    # Each iteration runs a short optimization, evaluates per-molecule
    # contributions, and randomly re-places the worst molecules.
    println("Packing $nfree free molecules ($(packmol_system.nmols) total)...")
    # Evaluate and print initial function value
    g0 = similar(x)
    return g0, x, cl_system, packmol_system, atom_positions, free_mol_indices
    f0 = fg!(g0, x, cl_system, packmol_system, atom_positions, free_mol_indices)
    @printf("Initial function value: %20.10e\n", f0)
    bestf = typemax(T)
    flast = typemax(T)
    converged = false
    best_positions = copy(packmol_system.molecule_positions)
    for loop in 0:nloop

        @printf("\n --- Starting packing loop: %d\n\n", loop)

        # Run a short optimization (maxit iterations per loop)
        spgresult = spgbox!(
            (g, x) -> fg!(g, x, cl_system, packmol_system, atom_positions, free_mol_indices),
            x;
            callback=(result) -> packmol_callback(result, cl_system, tol, iprint),
            vaux=auxvecs,
            nitmax=maxit,
            nfevalmax=10 * maxit,
        )

        # Update molecule positions with optimized values
        x_mol = reinterpret(MoleculePosition{D,T}, x)
        for (k, imol) in enumerate(free_mol_indices)
            packmol_system.molecule_positions[imol] = x_mol[k]
        end

        # Statistics: compute improvement of this loop relative to flast
        fx = spgresult.f
        dmin = min(cl_system.fg.dmin, cl_system.cutoff)
        fimp_last = flast > zero(T) ? clamp(-100 * (fx - flast) / flast, T(-99.99), T(99.99)) : T(100)
        fimprov = bestf < typemax(T) ? clamp(-100 * (fx - bestf) / bestf, T(-99.99), T(99.99)) : T(100)
        improved = fx < bestf
        if improved
            bestf = fx
            copyto!(best_positions, packmol_system.molecule_positions)
        end
        flast = fx

        @printf("\n  Loop %4d: f = %10.4e  best f = %10.4e", loop, fx, bestf)
        if loop > 0
            @printf("  improvement from best: %6.2f%%  from last: %6.2f%%", fimprov, fimp_last)
        end
        @printf("  min dist: %12.6f\n", dmin)

        # Check convergence: minimum distance exceeds tolerance
        if (tol - dmin < packmol_system.tolerance_precision) && bestf < packmol_system.constraint_precision
            print("""  
                Solution found at loop $(@sprintf("%d", loop))! 
                Minimum distance: $(@sprintf("%.6f", dmin))

            """)
            converged = true
            break
        end

        # Write best solution so far to output file
        if improved && !isempty(packmol_system.output_file)
            saved_positions = packmol_system.molecule_positions
            packmol_system.molecule_positions = best_positions
            write_output(packmol_system)
            packmol_system.molecule_positions = saved_positions
            println("  Best solution written to file: ", packmol_system.output_file)
        end

        # Move bad molecules if this loop did not improve f by at least 10%
        if fimp_last < T(10)
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
        @printf("  WARNING: packing did not converge after %d loops (best f = %.4e)\n", nloop, bestf)
    end

    # Restore best molecule positions
    copyto!(packmol_system.molecule_positions, best_positions)

    println("Minimum distance obtained: ", min(cl_system.fg.dmin, cl_system.cutoff))

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
    end

    tend = time()
    println("Packmol finished in ", (tend - tstart)," seconds." )

    return nothing
end

#
# SPGBox callback: print progress and check convergence
#
function packmol_callback(spgresult, cl_system, tol, iprint)
    if spgresult.nit % iprint == 0
        @printf(
            "  - Iteration: %6d  Minimum distance: %12.6f  Function value: %12.6e\n",
            spgresult.nit, min(cl_system.fg.dmin, cl_system.cutoff), spgresult.f
        )
    end
    if cl_system.fg.dmin > tol
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
    dmin = Inf
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
    @test dmin >= sys.tolerance - 0.1

    # Write output and verify
    output_file = Packmol.write_output(sys)
    @test isfile(output_file)
    output_atoms = read_pdb(output_file)
    @test length(output_atoms) == 300
    rm(output_file; force=true)
end