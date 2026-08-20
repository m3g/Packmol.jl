@testitem "_resnum" begin
    _resnum = Packmol._resnum
    # mode 0: sequential per molecule; fixed molecules always get 1
    @test _resnum(0, false, 1, 999, 1, 1, 1) == 1
    @test _resnum(0, false, 5, 999, 1, 1, 5) == 5
    @test _resnum(0, true, 1, 999, 1, 1, 1) == 1
    @test _resnum(0, true, 1, 999, 1, 1, 7) == 1
    # mode 0 wraps at 9999, with 0 mapped to 9999
    @test _resnum(0, false, 9999, 0, 1, 1, 1) == 9999
    @test _resnum(0, false, 10000, 0, 1, 1, 1) == 1
    # mode 1: verbatim template residue number, regardless of molecule index
    @test _resnum(1, false, 1, 42, 1, 1, 1) == 42
    @test _resnum(1, false, 5, 42, 1, 1, 5) == 42
    @test _resnum(1, false, 1, 0, 1, 1, 1) == 9999  # 0 still maps to 9999
    # mode 2: template offset continued by a running counter (irescount)
    @test _resnum(2, false, 1, 3, 1, 1, 1) == 3      # first molecule: no offset yet (irescount=1)
    @test _resnum(2, false, 2, 3, 1, 4, 2) == 6      # second molecule: irescount advanced by nres
    # mode 3: global molecule index across the whole system
    @test _resnum(3, false, 1, 999, 1, 1, 12) == 12
    @test _resnum(3, true, 1, 999, 1, 1, 12) == 12
    @test_throws ErrorException Packmol._resnum(4, false, 1, 1, 1, 1, 1)
end

@testitem "_chainc" begin
    _chainc = Packmol._chainc
    @test _chainc(0) == ' '
    @test _chainc(1) == 'A'
    @test _chainc(26) == 'Z'
    @test _chainc(27) == '1'
    @test _chainc(35) == '9'
    @test _chainc(36) == '0'
    @test _chainc(37) == '#'
    @test _chainc(1000) == '#'
end

@testitem "_unitcell_abc_angles" begin
    using Packmol: unitcell_matrix, _unitcell_abc_angles
    # Orthorhombic: angles must come out as 90/90/90
    uc = unitcell_matrix(Float64, 10.0, 20.0, 30.0, 90.0, 90.0, 90.0)
    a, b, c, α, β, γ = _unitcell_abc_angles(uc)
    @test a ≈ 10.0 && b ≈ 20.0 && c ≈ 30.0
    @test α ≈ 90.0 && β ≈ 90.0 && γ ≈ 90.0
    # Triclinic: round-trips through unitcell_matrix and back
    uc2 = unitcell_matrix(Float64, 12.0, 15.0, 18.0, 80.0, 95.0, 70.0)
    a2, b2, c2, α2, β2, γ2 = _unitcell_abc_angles(uc2)
    @test a2 ≈ 12.0 && b2 ≈ 15.0 && c2 ≈ 18.0
    @test α2 ≈ 80.0 && β2 ≈ 95.0 && γ2 ≈ 70.0
end

@testitem "write_output: residue numbering and chains" begin
    using PDBTools: read_pdb

    water_pdb = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")

    # Default (unset): single-residue template -> mode 0 (sequential per molecule)
    water = structure_type(water_pdb; number=5)
    sys = PackmolSystem([water]; output=tempname() * ".pdb", tolerance=2.0)
    outfile = write_output(sys)
    atoms = read_pdb(outfile)
    @test [a.resnum for a in atoms[1:3:end]] == collect(1:5)
    # Fewer than 9999 molecules of one type -> a single auto-assigned chain
    @test length(unique(a.chain for a in atoms)) == 1
    rm(outfile; force=true)

    # Explicit resnumbers=1: verbatim from the template (all molecules share
    # the template's own residue number, since water.pdb's is constant)
    water_verbatim = structure_type(water_pdb; number=3, residue_numbering=1)
    sys = PackmolSystem([water_verbatim]; output=tempname() * ".pdb", tolerance=2.0)
    outfile = write_output(sys)
    atoms = read_pdb(outfile)
    @test length(unique(a.resnum for a in atoms)) == 1
    rm(outfile; force=true)

    # Explicit resnumbers=3: global molecule index, continuing across
    # structure types (unlike mode 0, which restarts per type)
    water_a = structure_type(water_pdb; number=2, residue_numbering=3)
    water_b = structure_type(water_pdb; number=2, residue_numbering=3)
    sys = PackmolSystem([water_a, water_b]; output=tempname() * ".pdb", tolerance=2.0)
    outfile = write_output(sys)
    atoms = read_pdb(outfile)
    @test [a.resnum for a in atoms[1:3:end]] == collect(1:4)
    rm(outfile; force=true)

    # Explicit chain overrides auto-assignment
    water_chain = structure_type(water_pdb; number=3, chain='Z')
    sys = PackmolSystem([water_chain]; output=tempname() * ".pdb", tolerance=2.0)
    outfile = write_output(sys)
    atoms = read_pdb(outfile)
    @test all(a.chain == "Z" for a in atoms)
    rm(outfile; force=true)

    # changechains alternates between two chain letters, molecule by molecule
    water_alt = structure_type(water_pdb; number=4, changechains=true)
    sys = PackmolSystem([water_alt]; output=tempname() * ".pdb", tolerance=2.0)
    outfile = write_output(sys)
    atoms = read_pdb(outfile)
    mol_chains = [atoms[i].chain for i in 1:3:length(atoms)]
    @test mol_chains == ["A", "B", "A", "B"]
    rm(outfile; force=true)

    # chain + changechains together is rejected, matching Fortran
    @test_throws ErrorException structure_type(water_pdb; number=2, chain='A', changechains=true)
end

@testitem "write_output: CRYST1" begin
    using PDBTools: read_pdb

    water_pdb = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    water = structure_type(water_pdb; number=2, constraints=[InsideBox([0, 0, 0], [10, 10, 10])])

    # Neither add_box_sides nor PBC: no CRYST1 line
    sys = PackmolSystem([water]; output=tempname() * ".pdb", tolerance=2.0)
    outfile = write_output(sys)
    @test !any(occursin("CRYST1", line) for line in eachline(outfile))
    rm(outfile; force=true)

    # add_box_sides without PBC: CRYST1 derived from coordinate extrema
    sys = PackmolSystem([water]; output=tempname() * ".pdb", tolerance=2.0, add_box_sides=true)
    outfile = write_output(sys)
    lines = readlines(outfile)
    cryst1 = only(filter(l -> startswith(l, "CRYST1"), lines))
    @test occursin("90.00", cryst1)
    rm(outfile; force=true)

    # Explicit PBC: CRYST1 reflects the real unitcell, no CRYST1 without
    # explicit add_box_sides being required (PBC alone is enough)
    sys = PackmolSystem([water]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=Packmol.unitcell_matrix(Float64, 15.0, 15.0, 15.0, 90.0, 90.0, 90.0),
        unitcell_center=zeros(3),
    )
    outfile = write_output(sys)
    lines = readlines(outfile)
    cryst1 = only(filter(l -> startswith(l, "CRYST1"), lines))
    @test occursin("15.000", cryst1)
    rm(outfile; force=true)
end

@testitem "restart file round-trip" begin
    using StaticArrays

    positions = [
        Packmol.MoleculePosition(SVector(1.0, 2.0, 3.0), SVector(0.1, 0.2, 0.3)),
        Packmol.MoleculePosition(SVector(-4.5, 0.0, 9.25), SVector(-1.5, 3.0, 0.0)),
    ]
    file = tempname()
    Packmol._write_restart(file, positions)
    @test isfile(file)
    read_back = Packmol._read_restart(file, 2, Packmol.MoleculePosition{3,Float64})
    @test read_back == positions
    rm(file; force=true)

    # Too few lines for the requested molecule count
    short_file = tempname()
    Packmol._write_restart(short_file, positions[1:1])
    @test_throws ErrorException Packmol._read_restart(short_file, 2, Packmol.MoleculePosition{3,Float64})
    rm(short_file; force=true)
end

@testitem "PDB-based restart" begin
    using PDBTools: read_pdb, coor
    using StaticArrays
    using LinearAlgebra: norm

    @test Packmol._is_pdb_path("foo.pdb")
    @test Packmol._is_pdb_path("foo.PDB")
    @test Packmol._is_pdb_path("/some/dir/foo.pdb")
    @test !Packmol._is_pdb_path("foo.txt")
    @test !Packmol._is_pdb_path("foo")

    water_pdb = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    water = structure_type(water_pdb; number=3, constraints=[InsideBox([0, 0, 0], [28, 28, 28])])
    sys = PackmolSystem([water]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=Packmol.unitcell_matrix(Float64, 30.0, 30.0, 30.0, 90.0, 90.0, 90.0),
        unitcell_center=zeros(3),
    )
    @test packmol(sys; nloop=20, maxit=100, iprint=1000, seed=7)
    source_pdb = write_output(sys)

    # `_align_molecule` recovers each molecule's (cm, angles) well enough
    # that re-applying them to the reference template reproduces the
    # observed atom positions (up to the PDB format's own ~1e-3 precision).
    refcoords = water.reference_coordinates
    natoms = water.natoms
    atoms = read_pdb(source_pdb)
    for i in 1:water.number_of_molecules
        observed = SVector{3,Float64}[SVector{3,Float64}(coor(atoms[(i-1)*natoms+j])[1:3]) for j in 1:natoms]
        mp = Packmol._align_molecule(refcoords, observed)
        reconstructed = [Packmol.eulermat(mp.angles) * refcoords[j] + mp.cm for j in 1:natoms]
        @test maximum(norm(reconstructed[j] - observed[j]) for j in 1:natoms) < 1e-2
    end

    # Both the file (`.pdb`, dispatched by extension) and in-memory
    # `Vector{<:Atom}` sources go through the same alignment and must agree.
    segments = [(water.number_of_molecules, natoms, refcoords)]
    positions_from_file = Packmol._restart_positions(source_pdb, segments, Packmol.MoleculePosition{3,Float64})
    positions_from_atoms = Packmol._restart_positions(atoms, segments, Packmol.MoleculePosition{3,Float64})
    @test positions_from_file == positions_from_atoms

    # A source with the wrong atom count is a clear, immediate error.
    wrong_segments = [(water.number_of_molecules + 1, natoms, refcoords)]
    @test_throws ErrorException Packmol._restart_positions(atoms, wrong_segments, Packmol.MoleculePosition{3,Float64})

    rm(sys.output_file; force=true)
end
