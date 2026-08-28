@testitem "PackmolSystem from Julia code" begin
    using PDBTools: read_pdb
    file = Packmol.src_dir * "/../test/structure_files/water.pdb"

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
    @test _parse_options(String, "connect", "yes", ("yes" => true, "no" => false)) == true
    @test _parse_options(String, "connect", "no", ("yes" => true, "no" => false)) == false
    @test_throws ArgumentError _parse_options(String, "connect", "nop", ("yes" => true, "no" => false))
end

@testitem "read_packmol_input" begin
    using Packmol: read_packmol_input
    file = Packmol.src_dir * "/../test/input_files/water_box.inp"
    sys = read_packmol_input(file)
    @test sys.filetype == "pdb"
    @test sys.output_file == "water_box.pdb"
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

@testitem "radscale / discale keyword" begin
    using Packmol: read_packmol_input
    dir = dirname(Packmol.src_dir * "/../test/input_files/water_box.inp")
    original = read(Packmol.src_dir * "/../test/input_files/water_box.inp", String)

    mktempdir() do tmp
        cp(joinpath(dir, "water.pdb"), joinpath(tmp, "water.pdb"))

        radscale_file = joinpath(tmp, "radscale.inp")
        write(radscale_file, original * "\nradscale 1.5\n")
        sys = read_packmol_input(radscale_file)
        @test sys.radscale == 1.5

        discale_file = joinpath(tmp, "discale.inp")
        write(discale_file, original * "\ndiscale 1.3\n")
        sys = @test_logs (:warn, r"discale") read_packmol_input(discale_file)
        @test sys.radscale == 1.3
    end
end

@testitem "per-structure keyword does not leak into the global setting" begin
    using Packmol: read_packmol_input
    # `restart_from`/`restart_to` are valid both as a global keyword and as a
    # per-structure one (inside `structure ... end structure`); a
    # per-structure occurrence must only set that structure type's field,
    # not the whole-system `PackmolSystem` field of the same name —
    # regression test for a parser bug where the global-keyword dispatch
    # table was consulted unconditionally, before checking whether the
    # current line was inside a structure block.
    dir = dirname(Packmol.src_dir * "/../test/input_files/water_box.inp")
    original = read(Packmol.src_dir * "/../test/input_files/water_box.inp", String)

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
