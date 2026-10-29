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
    using PDBToolsBase: read_pdb

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
    using PDBToolsBase: read_pdb

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

@testitem "get_atoms / write_output: periodic_boundary_style" begin
    using StaticArrays: SVector
    using PDBToolsBase: read_pdb, position

    water_pdb = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    water = structure_type(water_pdb; number=1)

    d = 100.0
    unitcell = Packmol.dodecahedral_unitcell(Float64, d)
    center = SVector(0.0, 0.0, 0.0)
    sys = PackmolSystem([water]; output=tempname() * ".pdb", tolerance=2.0,
        unitcell=unitcell, unitcell_center=center,
    )
    @test sys.periodic_boundary_style == :triclinic

    # Place the molecule's CM at a point whose triclinic and dodecahedral
    # wraps disagree (see the round-trip test above for why: near a corner
    # of the parallelepiped domain that isn't the true nearest image), so the
    # two styles are actually distinguishable in what follows.
    cm = SVector(95.0, 95.0, 60.0)
    sys.molecule_positions[1] = Packmol.MoleculePosition(cm, zero(SVector{3,Float64}))
    cm_tri = Packmol.wrap_to_center(cm, unitcell, center)
    cm_dodeca = Packmol.triclinic_to_dodecahedral(cm, unitcell, center)
    @test cm_tri != cm_dodeca

    # Unwrapped Cartesian positions, independent of any PBC style — the
    # baseline `get_atoms` is expected to rigidly shift by (wrapped_cm - cm).
    natoms = length(sys.atoms)
    unwrapped = Vector{SVector{3,Float64}}(undef, natoms)
    Packmol.compute_atom_positions!(unwrapped, sys.molecule_positions, sys)

    atoms_tri = Packmol.get_atoms(sys)
    @test length(atoms_tri) == natoms
    for (a, u) in zip(atoms_tri, unwrapped)
        @test position(a) ≈ u + (cm_tri - cm) atol = 1e-4
    end

    Packmol.triclinic_to_dodecahedral(sys)
    @test sys.periodic_boundary_style == :dodecahedral
    atoms_dodeca = Packmol.get_atoms(sys)
    for (a, u) in zip(atoms_dodeca, unwrapped)
        @test position(a) ≈ u + (cm_dodeca - cm) atol = 1e-4
    end
    # the two styles actually produced different coordinates for this atom
    @test !(position(atoms_tri[1]) ≈ position(atoms_dodeca[1]))

    # write_output honors the same style (up to the PDB text format's own
    # 3-decimal-place precision)
    outfile = write_output(sys)
    written = read_pdb(outfile)
    for (a, w) in zip(atoms_dodeca, written)
        @test position(a) ≈ position(w) atol = 1e-3
    end
    rm(outfile; force=true)

    Packmol.dodecahedral_to_triclinic(sys)
    @test sys.periodic_boundary_style == :triclinic
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
    using PDBToolsBase: read_pdb, position
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
    @test packmol(sys; nloop=20, maxit=100, iprint=1000, seed=7).status == :packing_ready
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

@testitem "AMBER TER records" begin
    using PDBToolsBase: read_pdb
    water = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    # A copy of the water template with a TER record after its second atom
    template = tempname() * ".pdb"
    lines = readlines(water)
    iatom = findall(l -> startswith(l, "ATOM") || startswith(l, "HETATM"), lines)
    write(template, join(vcat(lines[1:iatom[2]], "TER", lines[iatom[2]+1:end]), "\n") * "\n")
    st = structure_type(template; number=3, constraints=[InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])])
    output = tempname() * ".pdb"
    nter(file) = count(l -> startswith(l, "TER"), readlines(file))
    # TER after each molecule
    sys = PackmolSystem([st]; output, tolerance=2.0, add_amber_ter=true)
    write_output(sys)
    @test nter(output) == 3
    @test length(read_pdb(output)) == 9
    out_lines = readlines(output)
    @test startswith(out_lines[findlast(l -> startswith(l, "ATOM") || startswith(l, "HETATM"), out_lines)+1], "TER")
    # TER records preserved from the template, and both options together
    sys = PackmolSystem([st]; output, tolerance=2.0, amber_ter_preserve=true)
    write_output(sys)
    @test nter(output) == 3
    out_lines = readlines(output)
    @test startswith(out_lines[findfirst(startswith("TER"), out_lines)-2], "ATOM") ||
          startswith(out_lines[findfirst(startswith("TER"), out_lines)-2], "HETATM")
    sys = PackmolSystem([st]; output, tolerance=2.0, add_amber_ter=true, amber_ter_preserve=true)
    write_output(sys)
    @test nter(output) == 6
    # No TER records by default
    sys = PackmolSystem([st]; output, tolerance=2.0)
    write_output(sys)
    @test nter(output) == 0
    rm(output; force=true)
    rm(template; force=true)
end

@testitem "segid" begin
    using PDBToolsBase: read_pdb
    water = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    constraints = [InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])]
    output = tempname() * ".pdb"
    sys = PackmolSystem([
            structure_type(water; number=2, constraints, segid="WAT"),
            structure_type(water; number=2, constraints),
        ]; output, tolerance=2.0)
    write_output(sys)
    atoms = read_pdb(output)
    @test all(at.segname == "WAT" for at in atoms[1:6])
    @test all(at.segname != "WAT" for at in atoms[7:12])
    @test_throws ArgumentError structure_type(water; number=2, segid="TOOLONG")
    rm(output; force=true)
end

@testitem "_FORCED output of a failed packing" begin
    water = joinpath(Packmol.src_dir, "..", "test", "structure_files", "water.pdb")
    output = tempname() * ".pdb"
    # Two molecules do not fit in this box (as Fortran's water_box_failed.inp)
    st = structure_type(water; number=2, constraints=[InsideBox([0.0, 0.0, 0.0], [3.0, 3.0, 1.0])])
    sys = PackmolSystem([st]; output, tolerance=2.0, nloop=1, max_iter=2)
    best_positions = copy(sys.molecule_positions)
    @test packmol(sys; iprint=1000).status == :failed_packing
    @test isfile(output)
    @test isfile(output * "_FORCED")
    rm(output; force=true)
    rm(output * "_FORCED"; force=true)
    # No _FORCED output when the packing succeeds
    st = structure_type(water; number=2, constraints=[InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])])
    sys = PackmolSystem([st]; output, tolerance=2.0)
    @test packmol(sys; iprint=1000).status == :packing_ready
    @test isfile(output)
    @test !isfile(output * "_FORCED")
    rm(output; force=true)
end

@testitem "CONECT records" begin
    using PDBToolsBase: read_pdb
    files = joinpath(Packmol.src_dir, "..", "test", "structure_files")
    ethanol = joinpath(files, "ethanol.pdb")
    water = joinpath(files, "water.pdb")
    constraints = [InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])]
    output = tempname() * ".pdb"
    nconect(file) = count(startswith("CONECT"), readlines(file))
    # The connectivity of the template is repeated for every copy of the molecule
    st_water = structure_type(water; number=2, constraints)
    st_ethanol = structure_type(ethanol; number=3, constraints)
    @test length(st_ethanol.connectivity) == 9
    @test st_water.connectivity == [[3], [3], [1, 2]]
    sys = PackmolSystem([st_water, st_ethanol]; output, tolerance=2.0)
    write_output(sys)
    @test nconect(output) == 33
    s = read_pdb(output)
    @test length(s.connectivity) == 33
    @test s.connectivity[1:6] == [[3], [3], [1, 2], [6], [6], [4, 5]]
    @test s.connectivity[7] == [8, 10, 11, 12]
    @test s.connectivity[16] == [17, 19, 20, 21]
    @test s.connectivity[33] == [27]
    # Also written with TER records
    sys = PackmolSystem([st_water, st_ethanol]; output, tolerance=2.0, add_amber_ter=true)
    write_output(sys)
    @test nconect(output) == 33
    @test read_pdb(output).connectivity == s.connectivity
    # Disabled for the whole system, or for a structure type
    sys = PackmolSystem([st_water, st_ethanol]; output, tolerance=2.0, connect=false)
    write_output(sys)
    @test nconect(output) == 0
    sys = PackmolSystem([st_water, structure_type(ethanol; number=3, constraints, connect=false)]; output, tolerance=2.0)
    write_output(sys)
    @test nconect(output) == 6
    rm(output; force=true)
end

@testitem "xyz files" begin
    using PDBToolsBase: read_xyz, read_pdb
    using Packmol: read_packmol_input
    files = joinpath(Packmol.src_dir, "..", "test", "structure_files")
    dir = mktempdir()
    input = joinpath(dir, "test.inp")
    # filetype xyz: the output is an xyz file, whatever its extension
    write(input, """
    tolerance 2.0
    filetype xyz
    output out.xyz
    structure $(joinpath(files, "ethanol.xyz"))
      number 3
      inside box 0. 0. 0. 20. 20. 20.
    end structure
    structure $(joinpath(files, "ethanol.pdb"))
      number 2
      connect no
      inside box 0. 0. 0. 20. 20. 20.
    end structure
    """)
    sys = read_packmol_input(input)
    @test sys.filetype == "xyz"
    @test sys.structure_types[1].natoms == 9
    @test isnothing(sys.structure_types[1].connectivity)
    @test isnothing(sys.structure_types[2].connectivity)
    write_output(sys)
    s = read_xyz(joinpath(dir, "out.xyz"))
    @test length(s) == 45
    @test s.title == "Built with Packmol"
    # filetype pdb: xyz structure files are recognized by the extension
    sys.filetype = "pdb"
    sys.output_file = joinpath(dir, "out.pdb")
    write_output(sys)
    @test length(read_pdb(joinpath(dir, "out.pdb"))) == 45
    # Unsupported file types
    write(input, "tolerance 2.0\nfiletype tinker\n")
    @test_throws ArgumentError read_packmol_input(input)
    @test_throws ArgumentError Packmol._file_format("water.txyz", "tinker")
    @test Packmol._file_format("water.txyz", "xyz") == "xyz"
    @test Packmol._file_format("water.XYZ", "pdb") == "xyz"
    @test Packmol._file_format("water.cif", "pdb") == "mmcif"
end
