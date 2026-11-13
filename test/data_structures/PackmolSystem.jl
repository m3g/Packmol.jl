@testitem "PackmolSystem from Julia code" setup=[TestFiles] begin
    using PDBToolsBase: read_pdb
    file = TestFiles.test_file("structure_files/water.pdb")

    st = structure_type(file; number=100, constraints=[InsideBox([0.,0.,0.],[40.,40.,40.])])
    sys = PackmolSystem([st]; output="julia_api_test.pdb", tolerance=2.0, seed=42)

    @test sys.output_file == "julia_api_test.pdb"
    @test sys.input_file == ""
    @test sys.tolerance == 2.0
    @test sys.seed == 42
    @test sys.nmols == 100
    @test length(sys.atoms) == 300
    @test all(at.molecule_index == 1 for at in sys.atoms[1:3])
    @test all(at.molecule_index == 100 for at in sys.atoms[298:300])
    @test length(sys.molecule_positions) == 100

    packmol(sys; nloop=200, maxit=20)
    @test isfile("julia_api_test.pdb")
    output_atoms = read_pdb("julia_api_test.pdb")
    @test length(output_atoms) == 300
    rm("julia_api_test.pdb"; force=true)
end

@testitem "fixed structure without number" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    water = TestFiles.test_file("structure_files/water.pdb")
    input = tempname() * ".inp"
    write(input, """
    tolerance 2.0
    structure $water
      number 5
      inside box 0. 0. 0. 20. 20. 20.
    end structure
    structure $water
      fixed 10. 10. 10. 0. 0. 0.
    end structure
    """)
    sys = read_packmol_input(input)
    @test sys.structure_types[2].number_of_molecules == 1
    @test sys.nmols == 6
    @test length(get_atoms(sys)) == 18
end

@testitem "fixed molecules: positions and periodic cell" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    using PDBToolsBase: position, read_pdb
    using StaticArrays: SVector, SMatrix
    using Statistics: mean
    # Fixed molecules are placed as in Fortran Packmol: rotated by the Euler angles
    # (`eulerfixed` of Fortran, `eulermat` here) about their geometric center with
    # `center`, or about the origin of their coordinates otherwise, and then translated.
    # The positions are the same in the output, in the packing, and in the overlap checks.
    water = TestFiles.test_file("structure_files/water.pdb")
    template = [SVector{3,Float64}(position(at)) for at in read_pdb(water)]
    cm, angles = SVector(5.0, -3.0, 2.0), SVector(0.7, 1.1, -0.4)
    R = Packmol.eulermat(angles)
    for (center, expected) in (
        (:geometric, [R * (x - mean(template)) + cm for x in template]),
        (:none, [R * x + cm for x in template]),
    )
        constraints = [InsideBox([-15.0, -15.0, -15.0], [15.0, 15.0, 15.0])]
        sys = PackmolSystem([
                structure_type(water; number=1, fixed=(cm, angles), center),
                structure_type(water; number=5, constraints),
            ]; output="", tolerance=2.0)
        redirect_stdout(devnull) do
            packmol(sys)
        end
        @test sys.status == :packing_ready
        @test all(isapprox(position(at), x; atol=1e-3) for (at, x) in zip(get_atoms(sys)[1:3], expected))
        @test all(isapprox(x, y; atol=1e-10) for (x, y) in zip(Packmol._collect_fixed_positions(sys), expected))
    end
    # Fixed molecules outside the periodic cell are an error (as Fortran Packmol's
    # testing/input_files/protein_outside_pbc_error.inp)
    @test_throws ArgumentError read_packmol_input(TestFiles.test_file("input_files/protein_outside_pbc_error.inp"))
    sys = read_packmol_input(TestFiles.test_file("input_files/solvprotein_pbc.inp"))
    @test isnothing(Packmol._check_fixed_inside_cell(sys))
    sys.structure_types[1].reference_coordinates .+= Ref(SVector(100.0, 0.0, 0.0))
    @test_throws ArgumentError packmol(sys)
    # Dodecahedral cells: the molecule may be outside the parallelepiped of the cell
    # vectors, but inside the dodecahedron
    ucell = SMatrix{3,3}(Packmol.dodecahedral_unitcell(Float64, 30.0))
    center = SVector(0.0, 0.0, 0.0)
    x = SVector(0.0, 0.0, 20.0)
    @test any(f -> abs(f) > 0.5, inv(ucell) * x)
    @test Packmol._inside_cell(x, ucell, inv(ucell), center)
    @test !Packmol._inside_cell(SVector(16.0, 0.0, 0.0), ucell, inv(ucell), center)
    # Atoms on the faces of a triclinic cell, slightly outside by the rounding of the
    # coordinates in PDB files (0.001 Å), are inside
    ucell = SMatrix{3,3}(Packmol.unitcell_matrix(Float64, 27.35, 27.35, 29.534, 90.0, 90.0, 120.0))
    face = -0.5 * ucell[:, 2] - 0.5 * ucell[:, 3]
    @test Packmol._inside_cell(face + SVector(0.0, -0.0005, -0.0005), ucell, inv(ucell), center)
    @test !Packmol._inside_cell(face + SVector(0.0, -0.1, 0.0), ucell, inv(ucell), center)
end

@testitem "show(::PackmolSystem)" setup=[TestFiles] begin
    using ShowMethodTesting
    water = TestFiles.test_file("structure_files/water.pdb")
    protein = TestFiles.test_file("structure_files/protein.pdb")
    sys = PackmolSystem([
            structure_type(protein; number=1, fixed=([0., 0., 0.], [0., 0., 0.])),
            structure_type(water; number=100, constraints=[InsideSphere([0., 0., 0.], 30.0)]),
        ]; output="show_test.pdb", tolerance=2.0)
    s = repr(MIME"text/plain"(), sys)
    @test parse_show(sys) ≈ """
    PackmolSystem{3,Float64}
        Number of types of structures: 2
        Total number of molecules: 101
        Structure types:
        1. "protein.pdb" - n_atoms = 4246 - n_mols = 1 - fixed - n_constraints = 0
        2. "water.pdb" - n_atoms = 3 - n_mols = 100 - n_constraints = 1"""
    @test !endswith(s, "\n")
end

@testitem "PackmolSystem: unitcell without unitcell_center" setup=[TestFiles] begin
    using StaticArrays: SVector
    water = TestFiles.test_file("structure_files/water.pdb")
    sys = PackmolSystem([structure_type(water; number=10)]; output=tempname() * ".pdb",
        unitcell=dodecahedral_unitcell(Float64, 15.0))
    @test sys.unitcell_center == SVector(0.0, 0.0, 0.0)
    @test redirect_stdout(devnull) do
        packmol(sys)
    end.status == :packing_ready
end
