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

    packmol(sys; nloop=200, maxit=20, iprint=200)
    @test isfile("julia_api_test.pdb")
    output_atoms = read_pdb("julia_api_test.pdb")
    @test length(output_atoms) == 300
    rm("julia_api_test.pdb"; force=true)
end

@testitem "_parse_value" begin
    using Packmol: _parse_value
    @test _parse_value(Int, "max_iter", "100") == 100
    @test _parse_value(Float64, "tolerance", "2.0") == 2.0
    @test _parse_value(Float32, "tolerance", "2.0") == 2.0f0
    @test_throws ArgumentError _parse_value(Int, "max_iter", "100.0")
end

@testitem "_parse_options" begin
    using Packmol: _parse_options
    @test _parse_options("connect", "yes", ("yes" => true, "no" => false)) == true
    @test _parse_options("connect", "no", ("yes" => true, "no" => false)) == false
    @test_throws ArgumentError _parse_options("connect", "nop", ("yes" => true, "no" => false))
end

@testitem "read_packmol_input" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    file = TestFiles.test_file("input_files/water_box.inp")
    sys = read_packmol_input(file)
    @test sys.filetype == "pdb"
    @test sys.output_file == joinpath(dirname(file), "water_box.pdb")
    @test sys.nmols == 1000
    @test length(sys.atoms) == 3000
    @test all(at.molecule_index == 1 for at in sys.atoms[1:3])
    @test all(at.molecule_index == 2 for at in sys.atoms[4:6])
    @test all(at.molecule_index == 1000 for at in sys.atoms[2998:3000])
    @test all(at.radius ≈ 1.0 for at in sys.atoms)
    @test length(sys.molecule_positions) == 1000
    @test sys.radscale == 1.2

    sys = read_packmol_input(file; T=Float32)
    @test typeof(sys.tolerance) == Float32
    @test eltype(sys.molecule_positions) == Packmol.MoleculePosition{3,Float32}
    @test sys.radscale ≈ 1.2

    # 2D: Currently not supported
    # sys = read_packmol_input(file; D=2)
end

@testitem "radscale / discale keyword" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    dir = dirname(TestFiles.test_file("input_files/water_box.inp"))
    original = read(TestFiles.test_file("input_files/water_box.inp"), String)

    mktempdir() do tmp
        cp(joinpath(dir, "water.pdb"), joinpath(tmp, "water.pdb"))

        radscale_file = joinpath(tmp, "radscale.inp")
        write(radscale_file, original * "\nradscale 1.5\n")
        sys = read_packmol_input(radscale_file)
        @test sys.radscale == 1.5

        discale_file = joinpath(tmp, "discale.inp")
        write(discale_file, original * "\ndiscale 1.3\n")
        sys = @test_warn r"discale" read_packmol_input(discale_file)
        @test sys.radscale == 1.3
    end
end

@testitem "per-structure keyword does not leak into the global setting" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    # `restart_from`/`restart_to` are valid both as a global keyword and as a
    # per-structure one (inside `structure ... end structure`); a
    # per-structure occurrence must only set that structure type's field,
    # not the whole-system `PackmolSystem` field of the same name —
    # regression test for a parser bug where the global-keyword dispatch
    # table was consulted unconditionally, before checking whether the
    # current line was inside a structure block.
    dir = dirname(TestFiles.test_file("input_files/water_box.inp"))
    original = read(TestFiles.test_file("input_files/water_box.inp"), String)

    mktempdir() do tmp
        cp(joinpath(dir, "water.pdb"), joinpath(tmp, "water.pdb"))
        restart_touched = joinpath(tmp, "some_restart_file.txt")

        text = replace(original, "end structure" => "restart_from $restart_touched\nend structure")
        file = joinpath(tmp, "per_structure_restart.inp")
        write(file, text)

        sys = read_packmol_input(file)
        @test isnothing(sys.restart_from)
        @test only(sys.structure_types).restart_from == restart_touched
    end
end

@testitem "read_packmol_input: comments, legacy keywords, nloop" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    water = TestFiles.test_file("structure_files/water.pdb")
    dir = mktempdir()
    input = joinpath(dir, "test.inp")
    write(input, """
    # A comment line
    tolerance 2.0 # inline comment
    output out.pdb
    nloop 37
    maxit 15
    movefrac 0.1
    movebadrandom
    short_tol_dist 1.0
    sidemax 1000.
    writecrd out.crd
    structure $water # the water
      number 10
      nloop 5
      segid WAT
      inside box 0. 0. 0. 20. 20. 20. # inline comment
    end structure
    """)
    sys = redirect_stderr(devnull) do
        read_packmol_input(input)
    end
    @test sys.tolerance == 2.0
    @test sys.nloop == 37
    @test sys.max_iter == 15
    @test sys.movefrac == 0.1
    @test sys.structure_types[1].segid == "WAT"
    @test sys.nmols == 10
    @test length(sys.structure_types[1].constraints) == 1
    # Not set in the input file: the `packmol` defaults apply
    sys = read_packmol_input(TestFiles.test_file("input_files/water_box.inp"))
    @test isnothing(sys.nloop)
    @test isnothing(sys.max_iter)
    @test isnothing(sys.movefrac)
    @test isnothing(sys.structure_types[1].segid)
    # segid has at most 4 characters
    write(input, "tolerance 2.0\nstructure $water\n  number 1\n  segid TOOLONG\nend structure\n")
    @test_throws ArgumentError read_packmol_input(input)
end

@testitem "read_packmol_input from IO" setup=[TestFiles] begin
    using Packmol: read_packmol_input
    # Paths are relative to the current directory
    dir = TestFiles.test_file("structure_files")
    input = "tolerance 2.0\noutput out.pdb\nstructure water.pdb\n  number 10\n  inside box 0. 0. 0. 20. 20. 20.\nend structure\n"
    sys = cd(dir) do
        read_packmol_input(IOBuffer(input))
    end
    @test sys.input_file == ""
    @test sys.output_file == "out.pdb"
    @test sys.nmols == 10
    @test sys.structure_types[1].filename == "water.pdb"
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
