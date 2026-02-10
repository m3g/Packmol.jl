using CellListMap: ParticleSystem, map_pairwise, map_pairwise!
using SPGBox: spgbox!

#
# Structure used to evaluate the function and gradient of the distance between atoms
# in a single pass, using CellListMap.
#
@kwdef mutable struct InteratomicDistanceFG{D,T}
    f::T
    g::Vector{MoleculePosition{D,T}}
    dmin::T
    fmol::Vector{T} # contribution of each molecule to the function value
    # Auxiliary array for gradient: carries the gradient relative to the cartesian coordinates of each atom
    gxcar::Vector{SVector{D,T}}
end
function InteratomicDistanceFG{D,T}(packmol_system::PackmolSystem) where {D,T}
    natoms = length(packmol_system.atoms)
    InteratomicDistanceFG(
        f = zero(T),
        g = zero(packmol_system.molecule_positions),
        dmin = typemax(T),
        fmol = zeros(T, packmol_system.nmols),
        gxcar = fill(zero(SVector{D,T}), natoms),
    )
end

# Custom copy, reset and reducer functions for CellListMap
function CellListMap.copy_output(x::InteratomicDistanceFG)
    InteratomicDistanceFG(
        x.f,
        copy(x.g),
        x.dmin,
        copy(x.fmol),
        copy(x.gxcar),
    )
end
function CellListMap.reset_output!(output::InteratomicDistanceFG{D,T}) where {D,T}
    output.f = zero(T)
    fill!(output.g, zero(MoleculePosition{D,T}))
    output.dmin = typemax(T)
    fill!(output.fmol, zero(T))
    fill!(output.gxcar, zero(SVector{D,T}))
    return output
end
function CellListMap.reducer(x::InteratomicDistanceFG, y::InteratomicDistanceFG)
    x.f += y.f
    x.g .+= y.g
    x.dmin = min(x.dmin, y.dmin)
    x.fmol .+= y.fmol
    x.gxcar .+= y.gxcar
    return x
end

# Updates the function and gradient of the system given a pair of
# particles within the cutoff.
function cartesian_fg!(x::T, y::T, i, j, d2, fg::InteratomicDistanceFG, packmol_system) where {T}
    iatom = packmol_system.atoms[i]
    jatom = packmol_system.atoms[j]
    if iatom.molecule_index == jatom.molecule_index
        return fg
    end
    tol = iatom.radius + jatom.radius
    d = sqrt(d2)
    fg.dmin = min(d, fg.dmin)
    if d < tol
        # Function value: add to total function value and to
        # molecular contributions
        fadd = (d - tol)^2
        fg.f += fadd
        fg.fmol[iatom.molecule_index] += fadd
        fg.fmol[jatom.molecule_index] += fadd
        # Gradient
        dv = y - x
        dvdd = 2 * (d - tol) * dv / d
        fg.gxcar[i] -= dvdd
        fg.gxcar[j] += dvdd
    end
    return fg
end

#
# Compute Cartesian atomic positions from rigid-body molecule DOFs.
# atom_position[i] = R(angles) * reference_coord[i] + center_of_mass
#
function compute_atom_positions!(
    atom_positions::Vector{SVector{D,T}},
    molecule_positions,
    packmol_system::PackmolSystem{D,T},
) where {D,T}
    iat = 0
    imol = 0
    for st in packmol_system.structure_types
        ref = st.reference_coordinates
        for _ in 1:st.number_of_molecules
            imol += 1
            mp = molecule_positions[imol]
            R = eulermat(mp.angles...)
            cm = mp.cm
            for j in 1:st.natoms
                iat += 1
                atom_positions[iat] = R * ref[j] + cm
            end
        end
    end
    return atom_positions
end

#
# Add constraint penalties and gradients to the fg output structure.
#
function constraint_fg!(
    fg_output::InteratomicDistanceFG{D,T},
    atom_positions::Vector{SVector{D,T}},
    packmol_system::PackmolSystem{D,T},
) where {D,T}
    for (iat, atom) in enumerate(packmol_system.atoms)
        x = atom_positions[iat]
        st = packmol_system.structure_types[atom.structure_type_index]
        for ic in atom.constraints
            c = st.constraints[ic]
            fg_output.f += constraint_penalty(c, x)
            fg_output.gxcar[iat] += constraint_gradient(c, x)
        end
    end
    return fg_output
end

#
# Combined objective function + gradient for SPGBox.
# Computes pairwise distance penalties + constraint penalties, then
# applies the chain rule to get gradients w.r.t. molecule DOFs.
#
function fg!(g, x,
    cl_system,
    packmol_system::PackmolSystem{D,T},
    atom_positions::Vector{SVector{D,T}},
) where {D,T}
    # Interpret optimizer variables as molecule positions and
    # update packmol_system so chain_rule! reads correct angles
    mol_pos = reinterpret(MoleculePosition{D,T}, x)
    packmol_system.molecule_positions .= mol_pos
    # Compute Cartesian atomic coordinates from molecule DOFs
    compute_atom_positions!(atom_positions, mol_pos, packmol_system)
    # Update CellListMap positions
    cl_system.xpositions .= atom_positions
    # Compute pairwise distance penalties and Cartesian gradients
    map_pairwise!(
        (x, y, i, j, d2, output) -> cartesian_fg!(x, y, i, j, d2, output, packmol_system),
        cl_system,
    )
    # Add constraint penalties and gradients
    constraint_fg!(cl_system.fg, atom_positions, packmol_system)
    # Chain rule: Cartesian → molecule DOF gradients
    chain_rule!(cl_system.fg, packmol_system)
    # Copy to SPGBox gradient vector
    g .= reinterpret(T, cl_system.fg.g)
    return cl_system.fg.f
end

#
# SPGBox callback: print progress and check convergence
#
function packmol_callback(spgresult, cl_system, tol, iprint)
    if spgresult.nit % iprint == 0
        @printf(
            " Iteration: %6d  Minimum distance: %12.6f  Function value: %12.6e\n",
            spgresult.nit, min(cl_system.fg.dmin, cl_system.cutoff), spgresult.f
        )
    end
    if cl_system.fg.dmin > tol
        return true
    end
    return false
end

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
    # Determine the placement region
    has_pbc = !isnothing(packmol_system.unitcell)
    if has_pbc
        # For PBC, place randomly within the unit cell
        uc = packmol_system.unitcell
    else
        # For non-PBC, place within the constraint bounding box
        lo, hi = compute_bounding_box(packmol_system)
    end
    imol = 0
    for st in packmol_system.structure_types
        for _ in 1:st.number_of_molecules
            imol += 1
            if st.fixed.fixed
                packmol_system.molecule_positions[imol] = st.fixed.position
            else
                if has_pbc
                    # Random fractional coordinates → Cartesian via unit cell matrix
                    frac = SVector{D,T}(ntuple(_ -> rand(RNG, T), D))
                    cm = SVector{D,T}(uc * frac)
                else
                    extent = hi - lo
                    cm = lo + SVector{D,T}(ntuple(d -> rand(RNG, T) * extent[d], D))
                end
                angles = SVector{D,T}(ntuple(_ -> T(2π) * rand(RNG, T), D))
                packmol_system.molecule_positions[imol] = MoleculePosition(cm, angles)
            end
        end
    end
    return packmol_system
end

"""
    packmol(packmol_system::PackmolSystem; kargs...)

Run the packing optimization on a `PackmolSystem`.
"""
function packmol(
    packmol_system::PackmolSystem{D,T};
    parallel::Bool=true,
    iprint::Int=10,
    nitmax::Int=1000,
    nfevalmax::Int=100*nitmax,
    seed::Int=packmol_system.seed,
) where {D,T}
    # Initialize RNG and molecule positions
    RNG = Random.Xoshiro(seed)
    println("Initializing molecule positions...")
    initialize_molecules!(packmol_system, RNG)

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
        lo, hi = compute_bounding_box(packmol_system)
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
    cutoff = max(packing_tol, (volume / ncells)^(one(T) / D))

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

    # Set up optimization variables (flat copy of molecule DOFs)
    x = copy(reinterpret(T, packmol_system.molecule_positions))
    auxvecs = SPGBox.VAux(x, zero(T))

    # Run optimization
    println("Packing $(packmol_system.nmols) molecules...")
    spgboxresult = spgbox!(
        (g, x) -> fg!(g, x, cl_system, packmol_system, atom_positions),
        x;
        callback=(spgresult) -> packmol_callback(spgresult, cl_system, tol, iprint),
        vaux=auxvecs,
        nitmax=nitmax,
        nfevalmax=nfevalmax,
    )

    # Update molecule positions with optimized values
    packmol_system.molecule_positions .= reinterpret(MoleculePosition{D,T}, x)

    println(spgboxresult)
    println("Minimum distance obtained: ", min(cl_system.fg.dmin, cl_system.cutoff))

    # For PBC: wrap atom positions into the first unit cell image
    if has_pbc
        compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)
        for i in eachindex(atom_positions)
            atom_positions[i] = CellListMap.wrap_to_first(atom_positions[i], cl_system.unitcell)
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

    return packmol_system
end

"""
    write_output(packmol_system::PackmolSystem)

Write packed coordinates to the output file specified in `packmol_system.output_file`.
"""
function write_output(packmol_system::PackmolSystem{D,T}) where {D,T}
    natoms = length(packmol_system.atoms)
    atom_positions = Vector{SVector{D,T}}(undef, natoms)
    compute_atom_positions!(atom_positions, packmol_system.molecule_positions, packmol_system)

    output_atoms = Atom[]
    iat = 0
    imol = 0
    atom_serial = 0
    for st in packmol_system.structure_types
        for _ in 1:st.number_of_molecules
            imol += 1
            for j in 1:st.natoms
                iat += 1
                atom_serial += 1
                atom = deepcopy(st.atoms[j])
                pos = atom_positions[iat]
                atom.index = atom_serial
                atom.x = pos[1]
                atom.y = pos[2]
                if D == 3
                    atom.z = pos[3]
                end
                atom.resnum = imol
                push!(output_atoms, atom)
            end
        end
    end

    output_file = packmol_system.output_file
    writePDB(output_atoms, output_file)
    println("Output written to: $output_file")
    return output_file
end

"""
    packmol(input_file::String; kargs...)

Read a Packmol input file, run the packing optimization, and write the output.
"""
function packmol(input_file::String; D::Int=3, T::DataType=Float64, kargs...)
    packmol_system = read_packmol_input(input_file; D, T)
    packmol(packmol_system; kargs...)
    return packmol_system
end

@testitem "monoatomic gradient" begin
    using StaticArrays
    using FiniteDifferences
    using CellListMap
    import Packmol: MonoAtomicFG, fg!

    # Testing function that computes the function value with the definition
    # of fg! above, to use finite-differences to check the gradient
    function f(x; dimension=2, unitcell=[1, 1], tol=0.1, parallel=false, return_grad=false)
        positions = [SVector{dimension}(x[i:i+dimension-1]) for i in 1:dimension:length(x)]
        system = ParticleSystem(
            xpositions=positions,
            unitcell=unitcell,
            cutoff=tol,
            output=Packmol.MonoAtomicFG(0.0, similar(positions), +Inf),
            output_name=:fg,
            parallel=parallel,
        )
        g = similar(x)
        if !return_grad
            return Packmol.fg!(g, x, system, tol)
        end
        f = Packmol.fg!(g, x, system, tol)
        return f, g
    end
    x = rand(2, 100)
    @test f(x; return_grad=true)[2] ≈ FiniteDifferences.grad(central_fdm(5, 1), f, x)[1] rtol = 1e-3
    @test f(x; return_grad=true, unitcell=[1 0.5; 0.5 1])[2] ≈
          FiniteDifferences.grad(central_fdm(5, 1), (x) ->
            f(x; unitcell=[1 0.5; 0.5 1]), x)[1] rtol = 1e-3
    x = rand(3, 100)
    @test f(x; dimension=3, unitcell=[1, 1, 1], return_grad=true)[2] ≈
          FiniteDifferences.grad(central_fdm(5, 1), (x) ->
            f(x; dimension=3, unitcell=[1, 1, 1]), x)[1] rtol = 1e-3
    @test f(x; dimension=3, unitcell=[1 0.5 0; 0 1 0.5; 0 0 1], return_grad=true)[2] ≈
          FiniteDifferences.grad(central_fdm(5, 1), (x) ->
            f(x; dimension=3, unitcell=[1 0.5 0; 0 1 0.5; 0 0 1]), x)[1] rtol = 1e-3

end # testitem gradient

@testitem "water box packing" begin
    using Packmol
    using PDBTools: readPDB
    using StaticArrays
    using LinearAlgebra: norm

    input_file = joinpath(Packmol.src_dir, "..", "test", "run_packmol", "water_box_small.inp")
    sys = Packmol.read_packmol_input(input_file)
    @test sys.nmols == 100
    @test length(sys.atoms) == 300
    @test sys.tolerance == 2.0

    Packmol.packmol(sys; nitmax=2000, iprint=50)

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
    output_atoms = readPDB(output_file)
    @test length(output_atoms) == 300
    rm(output_file; force=true)
end
