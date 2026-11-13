@testitem "water box packing" setup=[TestFiles] begin
    using Packmol
    using PDBToolsBase: read_pdb
    using StaticArrays
    using LinearAlgebra: norm

    input_file = TestFiles.test_file("input_files/water_box_small.inp")
    sys = Packmol.read_packmol_input(input_file)
    @test sys.nmols == 100
    @test length(sys.atoms) == 300
    @test sys.tolerance == 2.0
    @test sys.status == :not_packed

    result = Packmol.packmol(sys; nloop=200, maxit=20)
    @test result === sys
    @test sys.status == :packing_ready
    # Limited number of threads (restored afterwards)
    sys1 = Packmol.read_packmol_input(input_file)
    Packmol.packmol(sys1; nloop=200, maxit=20, nthreads=1)
    @test sys1.status == :packing_ready
    @test Packmol._ntasks[] == 0

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

@testitem "distance cutoff covers the radii in use" begin
    # Every penalized pair (d < r_i + r_j) must be within the cutoff, for the
    # working (radscale-inflated) radii and for the floor radii alike.
    floor = [1.0, 1.0, 1.5]
    @test Packmol._distance_cutoff(1.2 .* floor, floor, 2.0, 10.0, 3) ≈ 2 * 1.2 * 1.5
    @test Packmol._distance_cutoff(0.9 .* floor, floor, 2.0, 10.0, 3) ≈ 2 * 1.5
    # Never smaller than the tolerance, which `dmin` is compared with
    @test Packmol._distance_cutoff(0.5 .* floor, 0.5 .* floor, 2.0, 10.0, 3) ≈ 2.0
    # Grown when the box is so sparse that the cell count would exceed the cap
    @test Packmol._distance_cutoff(floor, floor, 2.0, 1e9, 3) ≈ (1e9 / (Packmol._MAX_CELLS_PER_ATOM * 3))^(1 / 3)
    # dmin is only known up to the cutoff
    @test Packmol._format_dmin(1.95, 2.4) == "1.950000"
    @test Packmol._format_dmin(2.4, 2.4) == "> 2.400000"
end

@testitem "structure radius larger than tolerance/2" setup=[TestFiles] begin
    using StaticArrays
    using LinearAlgebra: norm, Diagonal
    # A minority structure with a larger radius in a dense system: pairs
    # involving it must be kept at r_i + r_j > tolerance, which requires the
    # cell-list cutoff to follow the radii rather than the global tolerance.
    water_pdb = TestFiles.test_file("structure_files/water.pdb")
    L = 33.0
    water = structure_type(water_pdb; number=900)
    big = structure_type(water_pdb; number=30, radius=1.5)
    sys = PackmolSystem([water, big]; output=tempname() * ".pdb", tolerance=2.0, seed=1,
        unitcell=Matrix(Diagonal([L, L, L])), unitcell_center=SVector(0.0, 0.0, 0.0),
    )
    packmol(sys)
    @test sys.status == :packing_ready
    x = Vector{SVector{3,Float64}}(undef, length(sys.atoms))
    Packmol.compute_atom_positions!(x, sys.molecule_positions, sys)
    worst = let worst = Inf
        for i in eachindex(x), j in i+1:length(x)
            ai, aj = sys.atoms[i], sys.atoms[j]
            ai.molecule_index == aj.molecule_index && continue
            d = x[j] - x[i]
            d = d - L * round.(d / L)
            worst = min(worst, norm(d) - (ai.radius + aj.radius))
        end
        worst
    end
    @test worst >= -0.1
end

@testitem "the converged configuration is the solution" begin
    using StaticArrays
    using CellListMap: neighborlist
    # A converged loop may have a larger function value than the best earlier one (that
    # had fewer but larger overlaps): the converged configuration must be returned, not
    # the one of lowest function value. Gases in a zeolite (hexagonal cell, fixed
    # crystal, radii 1.5): seeds 4, 15, 16 and 20 returned overlaps of up to 0.14 Å.
    dir = mktempdir()
    src = joinpath(pkgdir(Packmol), "docs", "src", "tutorials", "files", "zeolite")
    for f in readdir(src)
        cp(joinpath(src, f), joinpath(dir, f))
    end
    inp = joinpath(dir, "zeolite.inp")
    input = read(inp, String)
    for seed in (4, 15, 16, 20)
        write(inp, "seed $seed\n" * input)
        sys = redirect_stdout(devnull) do
            packmol(inp; nthreads=1)
        end
        @test sys.status == :packing_ready
        x = similar(sys.atoms, SVector{3,Float64})
        Packmol.compute_atom_positions!(x, sys.molecule_positions, sys)
        worst = maximum(neighborlist(; positions=x, cutoff=3.0, unitcell=sys.unitcell); init=0.0) do (i, j, d)
            ai, aj = sys.atoms[i], sys.atoms[j]
            ai.molecule_index == aj.molecule_index ? 0.0 : ai.radius + aj.radius - d
        end
        @test worst < sys.tolerance_precision
    end
end

@testitem "non-periodic cell is not inflated by fixed molecules" setup=[TestFiles] begin
    using StaticArrays
    using LinearAlgebra: norm
    # A large fixed molecule surrounded by free molecules: the cell must cover
    # the region of the free molecules (plus their own extent) and the atoms of
    # the fixed molecule, but must not be inflated by the radius of the fixed one.
    protein_pdb = TestFiles.test_file("structure_files/protein.pdb")
    water_pdb = TestFiles.test_file("structure_files/water.pdb")
    protein = structure_type(protein_pdb; number=1, fixed=(zeros(3), zeros(3)), center=:geometric)
    water = structure_type(water_pdb; number=2)
    sys = PackmolSystem([protein, water]; output="", tolerance=2.0)
    Packmol._center_reference_coordinates!(sys)
    rprot = maximum(norm, sys.structure_types[1].reference_coordinates)
    rwat = maximum(norm, sys.structure_types[2].reference_coordinates)
    x = Vector{SVector{3,Float64}}(undef, length(sys.atoms))
    # Free molecules at opposite corners, beyond the fixed molecule
    L = rprot + 5
    sys.molecule_positions[2] = Packmol.MoleculePosition(SVector(-L, -L, -L), zero(SVector{3,Float64}))
    sys.molecule_positions[3] = Packmol.MoleculePosition(SVector(L, L, L), zero(SVector{3,Float64}))
    Packmol.compute_atom_positions!(x, sys.molecule_positions, sys)
    @test all(Packmol._nonperiodic_unitcell(sys, x) .≈ 1.2 * (2L + 2 * (rwat + 2.0)))
    # Free molecules close to each other: the cell covers the fixed atoms
    sys.molecule_positions[2] = Packmol.MoleculePosition(SVector(-1.0, 0.0, 0.0), zero(SVector{3,Float64}))
    sys.molecule_positions[3] = Packmol.MoleculePosition(SVector(1.0, 0.0, 0.0), zero(SVector{3,Float64}))
    Packmol.compute_atom_positions!(x, sys.molecule_positions, sys)
    natp = sys.structure_types[1].natoms
    lo = minimum(v -> v[1], x[1:natp]); hi = maximum(v -> v[1], x[1:natp])
    @test Packmol._nonperiodic_unitcell(sys, x)[1] ≈ 1.2 * (hi - lo + 2 * 2.0)
end

@testitem "packmol(input_file): paths relative to the input file" setup=[TestFiles] begin
    using Packmol
    # Run from another directory and check where the output lands.
    input_dir = mktempdir()
    for file in ("water_box_small.inp", "water.pdb")
        cp(TestFiles.test_file(joinpath("input_files", file)), joinpath(input_dir, file))
    end
    current_dir = pwd()
    system = cd(mktempdir()) do
        packmol(joinpath(input_dir, "water_box_small.inp"))
    end
    @test system.status == :packing_ready
    @test system.output_file == joinpath(input_dir, "water_box_small.pdb")
    @test isfile(system.output_file)
    @test pwd() == current_dir
end

@testitem "initial approximation that does not satisfy the constraints" setup=[TestFiles] begin
    water = TestFiles.test_file("structure_files/water.pdb")
    output = tempname() * ".pdb"
    # Incompatible constraints: inside a box, and outside a larger box around it
    input(extra) = IOBuffer("""
        tolerance 2.0
        output $output
        $extra
        structure $water
          number 10
          inside box 0. 0. 0. 10. 10. 10.
          outside box -1. -1. -1. 11. 11. 11.
        end structure
        """)
    sys = redirect_stdout(devnull) do
        packmol(read_packmol_input(input("")); nloop=5)
    end
    @test sys.status == :failed_initial_approximation
    @test isfile(output) # the initial approximation, for inspection
    rm(output)
    sys = read_packmol_input(input("continue_failed_initial_packing"))
    @test sys.continue_failed_initial_packing
    sys = redirect_stdout(devnull) do
        packmol(sys; nloop=5)
    end
    @test sys.status == :failed_packing
end

@testitem "report of the overlapping pairs" setup = [TestFiles] begin
    water = TestFiles.test_file("structure_files/water.pdb")
    # Too many molecules for the box, one loop: overlaps between the two types are reported
    a = structure_type(water; number = 40, constraints = [InsideBox([0., 0., 0.], [8., 8., 8.])])
    b = structure_type(water; number = 40, constraints = [InsideBox([0., 0., 0.], [8., 8., 8.])])
    sys = PackmolSystem([a, b]; tolerance = 2.0, output = tempname() * ".pdb", nloop = 1)
    out = tempname()
    open(out, "w") do io
        redirect_stdout(io) do
            packmol(sys)
        end
    end
    log = read(out, String)
    @test occursin("Overlapping pairs of atoms:", log)
    @test occursin("1: water.pdb - 2: water.pdb", log)
    # The counts are those of the overlaps beyond tolerance_precision, by pair of types
    fg = Packmol.InteratomicDistanceFG{3,Float64}(sys)
    @test size(fg.overlaps) == (2, 2)
    @test !fg.count_overlaps[]
end

@testitem "bad molecules: those that do not satisfy the convergence criterion" setup=[TestFiles] begin
    using Packmol: InteratomicDistanceFG, compute_atom_positions!, initialize_molecules!, fg!,
        _mark_unconverged!, MoleculePosition
    using StaticArrays: SVector
    using CellListMap: ParticleSystem
    using Random
    water = TestFiles.test_file("structure_files/water.pdb")
    st = structure_type(water; number=4, constraints=[InsideBox([0, 0, 0], [40, 40, 40])])
    sys = PackmolSystem([st]; output=tempname() * ".pdb", tolerance=2.0)
    initialize_molecules!(sys, Random.Xoshiro(1))
    angles = SVector(0.0, 0.0, 0.0)
    # 1 and 2 overlap by 0.2 Å (more than tolerance_precision, 0.01 Å), 3 is 0.3 Å outside
    # its box, 4 is fine
    sys.molecule_positions[1] = MoleculePosition(SVector(10.0, 10.0, 10.0), angles)
    sys.molecule_positions[2] = MoleculePosition(SVector(10.0, 10.0, 11.8), angles)
    sys.molecule_positions[4] = MoleculePosition(SVector(30.0, 30.0, 30.0), angles)
    ap = Vector{SVector{3,Float64}}(undef, length(sys.atoms))
    compute_atom_positions!(ap, sys.molecule_positions, sys)
    zmax = maximum(ap[i][3] for i in 1:3)
    sys.molecule_positions[3] = MoleculePosition(SVector(20.0, 20.0, 40.3 + 10.0 - zmax), angles)
    compute_atom_positions!(ap, sys.molecule_positions, sys)
    cl = ParticleSystem(xpositions=copy(ap), unitcell=[100.0, 100.0, 100.0], cutoff=2.2,
        output=InteratomicDistanceFG{3,Float64}(sys), output_name=:fg, parallel=false)
    free = collect(1:sys.nmols)
    x = copy(reinterpret(Float64, sys.molecule_positions))
    g = zeros(length(x))
    unconverged = zeros(Bool, sys.nmols)
    # recorded only when the flag is set
    fg!(g, x, cl, sys, ap, free)
    @test isempty(cl.fg.unconverged_mols)
    cl.fg.count_overlaps[] = true
    fg!(g, x, cl, sys, ap, free)
    cl.fg.count_overlaps[] = false
    _mark_unconverged!(unconverged, cl.fg)
    @test cl.fg.max_overlap >= sys.tolerance_precision
    @test cl.fg.max_constraint_penalty >= sys.constraint_precision
    @test unconverged == [true, true, true, false]
end

@testitem "radscale decay" setup=[TestFiles] begin
    using Packmol
    # exponential decay from radscale to 1, within 1e-3 of the initial excess at decay_loops
    @test Packmol._radscale_at(0, 3.0, 50.0) == 3.0
    @test 1 < Packmol._radscale_at(25, 3.0, 50.0) < 3
    @test Packmol._radscale_at(49, 3.0, 50.0) - 1 < 3e-3
    @test Packmol._radscale_at(50, 3.0, 50.0) == 1.0
    @test Packmol._radscale_at(10, 1.0, 50.0) == 1.0
    # the keyword argument and the system field
    input_file = TestFiles.test_file("input_files/water_box_small.inp")
    sys = Packmol.read_packmol_input(input_file)
    sys.radscale = 2.0
    sys.radscale_decay_loops = 5
    @test packmol(sys; nloop=200, seed=3).status == :packing_ready
    sys = Packmol.read_packmol_input(input_file)
    sys.radscale = 2.0
    @test packmol(sys; nloop=200, seed=3, radscale_decay_loops=5).status == :packing_ready
end

@testitem "writeout: minimum time between writes of the output" setup = [TestFiles] begin
    water = TestFiles.test_file("structure_files/water.pdb")
    function run(; kwargs...)
        a = structure_type(water; number = 300, constraints = [InsideBox([0., 0., 0.], [21., 21., 21.])])
        output = tempname() * ".pdb"
        sys = PackmolSystem([a]; tolerance = 2.0, output)
        out = tempname()
        open(out, "w") do io
            redirect_stdout(io) do
                packmol(sys; seed = 1, maxit = 5, kwargs...)
            end
        end
        return sys, output, count("Current solution written", read(out, String))
    end
    # every improved loop writes the output (maxit = 5: several loops)
    sys, output, nwrites = run(writeout = 0)
    @test sys.status == :packing_ready && isfile(output) && nwrites > 0
    # no intermediate writes; the final solution is written
    sys, output, nwrites = run(writeout = 1e6)
    @test sys.status == :packing_ready && isfile(output) && nwrites == 0
end

@testitem "global options as keyword arguments of packmol" setup = [TestFiles] begin
    using Packmol
    input_file = TestFiles.test_file("input_files/water_box_small.inp")
    sys = read_packmol_input(input_file)
    out = tempname() * ".pdb"
    redirect_stdout(devnull) do
        packmol(sys; tolerance = 2.2, maxit = 25, nloop = 100, seed = 5, output = out,
            writeout = 0, radscale_decay_loops = 7, add_box_sides = true, tolerance_precision = 0.05)
    end
    @test sys.status == :packing_ready
    @test sys.tolerance == 2.2 && all(a.radius == 1.1 for a in sys.atoms)
    @test all(st.radii == fill(1.1, st.natoms) for st in sys.structure_types)
    @test sys.max_iter == 25 && sys.nloop == 100 && sys.seed == 5 && sys.writeout == 0
    @test sys.radscale_decay_loops == 7 && sys.add_box_sides && sys.tolerance_precision == 0.05
    @test sys.output_file == out && isfile(out)
    # atoms with their own radius keep it
    sys = read_packmol_input(input_file)
    st = sys.structure_types[1]
    sys.atoms[1] = Packmol.AtomData{Float64}(sys.atoms[1].molecule_index, 1, 0.7, sys.atoms[1].constraints)
    Packmol._set_global_options!(sys; tolerance = 3.0)
    @test sys.atoms[1].radius == 0.7 && sys.atoms[2].radius == 1.5
    # periodic cells, as the input keywords
    Packmol._set_global_options!(sys; pbc = [30, 30, 30])
    @test sys.unitcell == [30 0 0; 0 30 0; 0 0 30] && sys.unitcell_center == [15, 15, 15]
    Packmol._set_global_options!(sys; pbc = (:octahedral, 0, 0, 0, 60))
    @test sys.wigner_seitz_cell && sys.unitcell_center == [0, 0, 0]
    Packmol._set_global_options!(sys; periodic_boundary_style = :octahedral)
    @test sys.periodic_boundary_style == :octahedral
    @test_throws ArgumentError Packmol._set_global_options!(sys; pbc = [1, 2])
end

@testitem "screen_log" setup=[TestFiles] begin
    using Packmol
    input_file = TestFiles.test_file("input_files/water_box.inp")
    # The report goes to stdout, with the progress bars appended (nothing redrawn)
    log = tempname()
    sys = open(log, "w") do io
        redirect_stdout(() -> packmol(input_file; nloop = 20), io)
    end
    report = read(log, String)
    @test sys.status == :packing_ready
    @test occursin("|*****", report)
    @test !occursin('\r', report)
    # `screen_log = false`: nothing printed
    sys = open(log, "w") do io
        redirect_stdout(() -> packmol(input_file; nloop = 20, screen_log = false), io)
    end
    @test sys.status == :packing_ready
    @test isempty(read(log, String))
    @test Packmol._screen_log[] # restored after the run
    # The input keyword, which also silences the warnings of the parser
    quiet_file = joinpath(dirname(input_file), "water_box_quiet.inp")
    write(quiet_file, "screen_log no\nfscale 1.0\n" * read(input_file, String))
    err = tempname()
    sys = open(err, "w") do io
        redirect_stderr(() -> read_packmol_input(quiet_file), io)
    end
    @test !sys.screen_log
    @test isempty(read(err, String))
end

@testitem "same seed, one task: the same packing" setup=[TestFiles] begin
    # The partial sums of the threaded loops are added in a fixed order (not as the tasks
    # finish), so with one task (and the constraint loop still split into several tasks)
    # two runs give the same result. With several tasks, the pairs of CellListMap are
    # distributed dynamically, and the results can differ in the last digits.
    using Packmol
    input_file = TestFiles.test_file("input_files/bilayer_pbc.inp")
    outputs = map(1:2) do run
        output = joinpath(mktempdir(), "bilayer_pbc.pdb")
        redirect_stdout(devnull) do
            packmol(input_file; output, nthreads = 1)
        end
        read(output, String)
    end
    @test outputs[1] == outputs[2]
end
