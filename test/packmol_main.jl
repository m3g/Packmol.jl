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

    result = Packmol.packmol(sys; nloop=200, maxit=20, iprint=10)
    @test result === sys
    @test sys.status == :packing_ready
    # Limited number of threads (restored afterwards)
    sys1 = Packmol.read_packmol_input(input_file)
    Packmol.packmol(sys1; nloop=200, maxit=20, iprint=1000, nthreads=1)
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

@testitem "constrain_rotation" setup=[TestFiles] begin
    using LinearAlgebra: norm, dot
    using StaticArrays: SVector
    using PDBToolsBase: position
    angle_z(v) = acosd(clamp(v[3] / norm(v), -1, 1))
    lipid = TestFiles.test_file("structure_files/palmitoil.pdb")
    # The vector from the center of atoms 1 and 2 to atom 42, within 20° of z
    st = structure_type(lipid; number=20,
        constraints=[InsideBox([-20, -20, -20], [20, 20, 20])],
        constrain_rotation=Dict(:z => (0.0, 20.0, [1, 2], [42])),
    )
    sys = PackmolSystem([st]; output=tempname() * ".pdb", tolerance=2.0)
    packmol(sys; iprint=1000)
    @test sys.status == :packing_ready
    x = [SVector{3,Float64}(position(at)) for at in get_atoms(sys)]
    angles = [angle_z(x[(i-1)*42+42] - (x[(i-1)*42+1] + x[(i-1)*42+2]) / 2) for i in 1:20]
    @test maximum(angles) < 21.0
    # Default vector: the longest axis of the molecule, here around 90° from x
    st = structure_type(lipid; number=20,
        constraints=[InsideBox([-20, -20, -20], [20, 20, 20])],
        constrain_rotation=Dict(:x => (90.0, 10.0)),
    )
    sys = PackmolSystem([st]; output=tempname() * ".pdb", tolerance=2.0)
    packmol(sys; iprint=1000)
    @test sys.status == :packing_ready
    rc = only(st.rotation_constraints)
    for mp in sys.molecule_positions
        v = Packmol.eulermat(mp.angles) * rc.reference_vector
        @test 79.0 < acosd(clamp(v[1] / norm(v), -1, 1)) < 101.0
    end
end

@testitem "rotation constraints: penalty, gradient and weights" setup=[TestFiles] begin
    using StaticArrays: SVector
    using LinearAlgebra: norm
    RC = Packmol.RotationConstraint{3,Float64}
    rc = RC(:z, 0.0, 10.0; from=[1], to=[2])
    @test rc.lower == 0.0 && rc.upper ≈ deg2rad(10)
    @test Packmol.rotation_penalty(rc, SVector(0.0, 0.1, 1.0)) == 0.0
    v = SVector(0.3, -0.5, 0.8)
    @test Packmol.rotation_penalty(rc, v) > 0
    # Gradient against finite differences
    h = 1e-6
    g = Packmol.rotation_gradient(rc, v)
    gnum = SVector(ntuple(i -> (Packmol.rotation_penalty(rc, v + h * SVector(ntuple(j -> j == i, 3))) -
        Packmol.rotation_penalty(rc, v - h * SVector(ntuple(j -> j == i, 3)))) / 2h, 3))
    # (at constant |v|, as for a rigid molecule: the component along v is not used)
    vhat = v / norm(v)
    @test g ≈ gnum - (gnum' * vhat) * vhat rtol = 1e-5
    # Lower bound (angle with x between 80° and 100°)
    rc = RC(:x, 90.0, 10.0)
    @test Packmol.rotation_penalty(rc, SVector(1.0, 0.0, 0.2)) > 0
    @test Packmol.rotation_penalty(rc, SVector(0.1, 0.0, 1.0)) == 0
    # Default weights: v = Σ wᵢ rᵢ is the z axis of the reference orientation, with the
    # length of the molecule along it
    ref = [SVector(cos(i), 0.5 * sin(i), 3.0 * i) for i in 1:10]
    Packmol.set_reference_coordinates!(ref)
    Packmol._set_rotation_constraint!(rc, ref)
    Lz = maximum(r -> r[3], ref) - minimum(r -> r[3], ref)
    @test rc.reference_vector ≈ SVector(0.0, 0.0, Lz)
    @test sum(rc.weights) ≈ 0 atol = 1e-12
    # Atom sets: the vector between their centers
    rc = RC(:z, 0.0, 10.0; from=[1, 2], to=[10])
    Packmol._set_rotation_constraint!(rc, ref)
    @test rc.reference_vector ≈ ref[10] - (ref[1] + ref[2]) / 2
    @test_throws ArgumentError RC(:w, 0.0, 10.0)
    @test_throws ArgumentError RC(:z, 0.0, 10.0; from=[1])
end

@testitem "rotation constraints: gradient of the objective function" setup=[TestFiles] begin
    using Packmol: InteratomicDistanceFG, compute_atom_positions!, initialize_molecules!, fg!
    using StaticArrays: SVector
    using CellListMap: ParticleSystem
    using Random
    lipid = TestFiles.test_file("structure_files/palmitoil.pdb")
    st = structure_type(lipid; number=4,
        constraints=[InsideBox([-10, -10, -10], [10, 10, 10])],
        constrain_rotation=Dict(:z => (0.0, 5.0, [1, 2], [40, 42]), :x => (90.0, 5.0)),
    )
    sys = PackmolSystem([st]; output=tempname() * ".pdb", tolerance=2.0)
    initialize_molecules!(sys, Random.Xoshiro(1))
    ap = Vector{SVector{3,Float64}}(undef, length(sys.atoms))
    compute_atom_positions!(ap, sys.molecule_positions, sys)
    cl = ParticleSystem(xpositions=copy(ap), unitcell=[100.0, 100.0, 100.0], cutoff=2.2,
        output=InteratomicDistanceFG{3,Float64}(sys), output_name=:fg, parallel=false)
    free = collect(1:sys.nmols)
    x = copy(reinterpret(Float64, sys.molecule_positions))
    g = zeros(length(x))
    fg!(g, x, cl, sys, ap, free)
    @test sys.structure_types[1].rotation_constraints[1].reference_vector != zero(SVector{3,Float64})
    h = 1e-6
    gnum = map(eachindex(x)) do i
        xp = copy(x); xp[i] += h
        xm = copy(x); xm[i] -= h
        (fg!(zeros(length(x)), xp, cl, sys, ap, free) - fg!(zeros(length(x)), xm, cl, sys, ap, free)) / 2h
    end
    @test g ≈ gnum rtol = 1e-4
end

@testitem "restart_from / restart_to" setup=[TestFiles] begin
    water_pdb = TestFiles.test_file("structure_files/water.pdb")
    urea_pdb = TestFiles.test_file("structure_files/urea.pdb")
    uc = Packmol.unitcell_matrix(Float64, 30.0, 30.0, 30.0, 90.0, 90.0, 90.0)

    restart_all = tempname() * ".pdb"
    restart_water = tempname() * ".xyz"

    # First run: converge normally, writing the coordinates of the whole system (PDB)
    # and of the water molecules (XYZ)
    water = structure_type(water_pdb; number=4,
        constraints=[InsideBox([0, 0, 0], [28, 28, 28])], restart_to=restart_water)
    urea = structure_type(urea_pdb; number=2, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys1 = PackmolSystem([water, urea]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3), restart_to=restart_all,
    )
    @test packmol(sys1; nloop=20, maxit=100, iprint=1000, seed=1).status == :packing_ready
    @test isfile(restart_all)
    @test isfile(restart_water)
    @test countlines(restart_water) == 2 + 4 * 3
    saved_water_cm = sys1.molecule_positions[1].cm
    rm(sys1.output_file; force=true)

    # Second run: whole-system restart_from should skip initial placement
    # entirely and reproduce the exact same (already-converged) positions —
    # converging immediately, in the very first loop.
    water2 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    urea2 = structure_type(urea_pdb; number=2, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys2 = PackmolSystem([water2, urea2]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3), restart_from=restart_all,
    )
    @test packmol(sys2; nloop=20, maxit=100, iprint=1000, seed=2).status == :packing_ready
    @test sys2.molecule_positions[1].cm ≈ saved_water_cm atol = 1e-2
    rm(sys2.output_file; force=true)

    # Third run: per-structure-type restart_from applies only to that type,
    # and (regression test) must NOT leak into the whole-system setting —
    # the other structure type still gets a normal (non-restarted) placement.
    water3 = structure_type(water_pdb; number=4,
        constraints=[InsideBox([0, 0, 0], [28, 28, 28])], restart_from=restart_water)
    urea3 = structure_type(urea_pdb; number=2, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys3 = PackmolSystem([water3, urea3]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3),
    )
    @test isnothing(sys3.restart_from)
    @test sys3.structure_types[1].restart_from == restart_water
    @test packmol(sys3; nloop=20, maxit=100, iprint=1000, seed=3).status == :packing_ready
    @test sys3.molecule_positions[1].cm ≈ saved_water_cm atol = 1e-2
    rm(sys3.output_file; force=true)

    rm(restart_all; force=true)
    rm(restart_water; force=true)
end

@testitem "restart_from a PDB file" setup=[TestFiles] begin
    using PDBToolsBase: read_pdb

    water_pdb = TestFiles.test_file("structure_files/water.pdb")
    uc = Packmol.unitcell_matrix(Float64, 30.0, 30.0, 30.0, 90.0, 90.0, 90.0)

    # First run: a normal PDB output file, no restart involved — it's not
    # Packmol's own restart format, but it should work as a restart_from
    # source anyway, since it's just atom positions like any other PDB.
    water1 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys1 = PackmolSystem([water1]; output=tempname() * ".pdb", tolerance=2.0, unitcell=uc, unitcell_center=zeros(3))
    @test packmol(sys1; nloop=20, maxit=100, iprint=1000, seed=11).status == :packing_ready
    source_pdb = sys1.output_file

    # Restarting from that PDB (file path, dispatched by its .pdb extension)
    # should reproduce essentially the same, already-converged positions —
    # converging immediately — using a different seed to rule out the run
    # just coincidentally reconverging to the same configuration on its own.
    water2 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys2 = PackmolSystem([water2]; output=tempname() * ".pdb", tolerance=2.0, unitcell=uc, unitcell_center=zeros(3),
        restart_from=source_pdb,
    )
    @test packmol(sys2; nloop=20, maxit=100, iprint=1000, seed=22).status == :packing_ready
    for i in 1:4
        @test sys2.molecule_positions[i].cm ≈ sys1.molecule_positions[i].cm atol = 1e-2
    end
    rm(sys2.output_file; force=true)

    # Same, but passing an already-loaded `Vector{<:Atom}` directly (the
    # Julia-API-only form) instead of a file path.
    atoms = read_pdb(source_pdb)
    water3 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys3 = PackmolSystem([water3]; output=tempname() * ".pdb", tolerance=2.0, unitcell=uc, unitcell_center=zeros(3),
        restart_from=atoms,
    )
    @test packmol(sys3; nloop=20, maxit=100, iprint=1000, seed=33).status == :packing_ready
    for i in 1:4
        @test sys3.molecule_positions[i].cm ≈ sys1.molecule_positions[i].cm atol = 1e-2
    end
    rm(sys3.output_file; force=true)

    rm(source_pdb; force=true)
end

@testitem "restart_from with a fixed molecule present" setup=[TestFiles] begin
    water_pdb = TestFiles.test_file("structure_files/water.pdb")
    uc = Packmol.unitcell_matrix(Float64, 30.0, 30.0, 30.0, 90.0, 90.0, 90.0)

    restart_all = tempname() * ".pdb"

    # Fixed structure type declared *before* the free one — matching common
    # real-world inputs (e.g. a fixed protein listed first, solvent after) —
    # so a whole-system restart must correctly skip over the fixed
    # molecule's own atoms/lines (interleaved at the front of the restart
    # source) instead of assuming the source only ever lists free molecules.
    fixed_water = structure_type(water_pdb; number=1, fixed=(zeros(3), zeros(3)))
    water1 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    # `avoid_overlap=false`: with only 3 fixed atoms, the fixed-particle
    # overlap-check system's own auto-sized cutoff (see
    # `_build_overlap_check_system`) can exceed half its own bounding box —
    # a pre-existing limitation for tiny fixed structures, unrelated to what
    # this test is actually checking (restart_from's atom bookkeeping).
    sys1 = PackmolSystem([fixed_water, water1]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3), restart_to=restart_all, avoid_overlap=false,
    )
    @test packmol(sys1; nloop=20, maxit=100, iprint=1000, seed=41).status == :packing_ready
    source_pdb = sys1.output_file
    saved_free_cm = [sys1.molecule_positions[i].cm for i in 2:5]

    # Whole-system restart_from a PDB written for the entire system (fixed
    # molecule included) must not error on the atom count, and must recover
    # the free molecules' positions (the fixed one keeps its own `fixed`
    # position regardless of what the restart source contains for it).
    fixed_water2 = structure_type(water_pdb; number=1, fixed=(zeros(3), zeros(3)))
    water2 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys2 = PackmolSystem([fixed_water2, water2]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3), restart_from=source_pdb, avoid_overlap=false,
    )
    @test packmol(sys2; nloop=20, maxit=100, iprint=1000, seed=42).status == :packing_ready
    for i in 2:5
        @test sys2.molecule_positions[i].cm ≈ saved_free_cm[i-1] atol = 1e-2
    end
    rm(sys2.output_file; force=true)

    # Same, but restarting from the whole-system restart_to file (also
    # written for the entire system, fixed molecule included).
    fixed_water3 = structure_type(water_pdb; number=1, fixed=(zeros(3), zeros(3)))
    water3 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys3 = PackmolSystem([fixed_water3, water3]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3), restart_from=restart_all, avoid_overlap=false,
    )
    @test packmol(sys3; nloop=20, maxit=100, iprint=1000, seed=43).status == :packing_ready
    for i in 2:5
        @test sys3.molecule_positions[i].cm ≈ saved_free_cm[i-1] atol = 1e-2
    end
    rm(sys3.output_file; force=true)

    rm(source_pdb; force=true)
    rm(restart_all; force=true)
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
        packmol(joinpath(input_dir, "water_box_small.inp"); iprint=1000)
    end
    @test system.status == :packing_ready
    @test system.output_file == joinpath(input_dir, "water_box_small.pdb")
    @test isfile(system.output_file)
    @test pwd() == current_dir
end
