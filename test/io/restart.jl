@testitem "restart_from" setup=[TestFiles] begin
    water_pdb = TestFiles.test_file("structure_files/water.pdb")
    urea_pdb = TestFiles.test_file("structure_files/urea.pdb")
    uc = Packmol.unitcell_matrix(Float64, 30.0, 30.0, 30.0, 90.0, 90.0, 90.0)

    restart_water = tempname() * ".xyz"

    # First run: converge normally; its output (PDB) is the restart file of the whole
    # system, and the coordinates of its water molecules (XYZ) that of the water
    water = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    urea = structure_type(urea_pdb; number=2, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys1 = PackmolSystem([water, urea]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3),
    )
    @test packmol(sys1; nloop=20, maxit=100, seed=1).status == :packing_ready
    restart_all = sys1.output_file
    Packmol.write_xyz(restart_water, Packmol.Structure(get_atoms(sys1)[1:12]))
    @test countlines(restart_water) == 2 + 4 * 3
    saved_water_cm = sys1.molecule_positions[1].cm

    # Second run: whole-system restart_from should skip initial placement
    # entirely and reproduce the exact same (already-converged) positions —
    # converging immediately, in the very first loop.
    water2 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    urea2 = structure_type(urea_pdb; number=2, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys2 = PackmolSystem([water2, urea2]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3), restart_from=restart_all,
    )
    @test packmol(sys2; nloop=20, maxit=100, seed=2).status == :packing_ready
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
    @test packmol(sys3; nloop=20, maxit=100, seed=3).status == :packing_ready
    # The urea is placed anew and packed around the restarted water, which can move a little
    # (0.24 Å on CI, depending on the platform and the number of threads): a water placed
    # anew would be anywhere in the box
    @test sys3.molecule_positions[1].cm ≈ saved_water_cm atol = 1.0
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
    @test packmol(sys1; nloop=20, maxit=100, seed=11).status == :packing_ready
    source_pdb = sys1.output_file

    # Restarting from that PDB (file path, dispatched by its .pdb extension)
    # should reproduce essentially the same, already-converged positions —
    # converging immediately — using a different seed to rule out the run
    # just coincidentally reconverging to the same configuration on its own.
    water2 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys2 = PackmolSystem([water2]; output=tempname() * ".pdb", tolerance=2.0, unitcell=uc, unitcell_center=zeros(3),
        restart_from=source_pdb,
    )
    @test packmol(sys2; nloop=20, maxit=100, seed=22).status == :packing_ready
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
    @test packmol(sys3; nloop=20, maxit=100, seed=33).status == :packing_ready
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
        unitcell=uc, unitcell_center=zeros(3), avoid_overlap=false,
    )
    @test packmol(sys1; nloop=20, maxit=100, seed=41).status == :packing_ready
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
    @test packmol(sys2; nloop=20, maxit=100, seed=42).status == :packing_ready
    for i in 2:5
        @test sys2.molecule_positions[i].cm ≈ saved_free_cm[i-1] atol = 1e-2
    end
    rm(sys2.output_file; force=true)

    # Same, but restarting from the output written as mmCIF (also for the
    # entire system, fixed molecule included).
    restart_all = tempname() * ".cif"
    write_output(sys1; output_file=restart_all)
    fixed_water3 = structure_type(water_pdb; number=1, fixed=(zeros(3), zeros(3)))
    water3 = structure_type(water_pdb; number=4, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys3 = PackmolSystem([fixed_water3, water3]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=uc, unitcell_center=zeros(3), restart_from=restart_all, avoid_overlap=false,
    )
    @test packmol(sys3; nloop=20, maxit=100, seed=43).status == :packing_ready
    for i in 2:5
        @test sys3.molecule_positions[i].cm ≈ saved_free_cm[i-1] atol = 1e-2
    end
    rm(sys3.output_file; force=true)

    rm(source_pdb; force=true)
    rm(restart_all; force=true)
end

@testitem "restart files are coordinate files" begin
    for file in ("a.pdb", "a.PDB", "dir/a.xyz", "a.cif", "a.mmcif")
        @test Packmol._check_restart_file(file, "restart_from") == file
    end
    # The restart files of Fortran Packmol are not supported
    @test_throws ArgumentError Packmol._check_restart_file("water.pack", "restart_from")
    @test_throws ArgumentError Packmol._check_restart_file("restart", "restart_from")
    @test isnothing(Packmol._check_restart_file(nothing, "restart_from"))
end

@testitem "PDB-based restart" setup=[TestFiles] begin
    using PDBToolsBase: read_pdb, position
    using StaticArrays
    using LinearAlgebra: norm

    water_pdb = TestFiles.test_file("structure_files/water.pdb")
    water = structure_type(water_pdb; number=3, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys = PackmolSystem([water]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=Packmol.unitcell_matrix(Float64, 30.0, 30.0, 30.0, 90.0, 90.0, 90.0),
        unitcell_center=zeros(3),
    )
    @test packmol(sys; nloop=20, maxit=100, seed=7).status == :packing_ready
    source_pdb = write_output(sys)

    # `_align_molecule` recovers each molecule's (cm, angles) well enough
    # that re-applying them to the reference template reproduces the
    # observed atom positions (up to the PDB format's own ~1e-3 precision).
    refcoords = water.reference_coordinates
    natoms = water.natoms
    atoms = read_pdb(source_pdb)
    for i in 1:water.number_of_molecules
        observed = SVector{3,Float64}[SVector{3,Float64}(position(atoms[(i-1)*natoms+j])[1:3]) for j in 1:natoms]
        mp = Packmol._align_molecule(refcoords, observed)
        reconstructed = [Packmol.eulermat(mp.angles) * refcoords[j] + mp.cm for j in 1:natoms]
        @test maximum(norm(reconstructed[j] - observed[j]) for j in 1:natoms) < 1e-2
    end

    # Both the file (`.pdb`, dispatched by extension) and in-memory
    # `Vector{<:Atom}` sources go through the same alignment and must agree.
    segments = [(water.number_of_molecules, natoms, refcoords, true)]
    positions_from_file = Packmol._restart_positions(source_pdb, segments, Packmol.MoleculePosition{3,Float64})
    positions_from_atoms = Packmol._restart_positions(atoms, segments, Packmol.MoleculePosition{3,Float64})
    @test positions_from_file == positions_from_atoms

    # A source with the wrong atom count is a clear, immediate error.
    wrong_segments = [(water.number_of_molecules + 1, natoms, refcoords, true)]
    @test_throws ErrorException Packmol._restart_positions(atoms, wrong_segments, Packmol.MoleculePosition{3,Float64})

    rm(sys.output_file; force=true)
end
