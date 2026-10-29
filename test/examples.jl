#
# Regression tests based on the classic Fortran Packmol example inputs
# (full molecule counts, not the scaled-down variants used elsewhere in
# test/input_files). Each one must actually converge, and in a number of
# packing loops not much greater than what it currently needs — a large
# increase would signal a regression in the stall-detection/movebad!
# machinery in packmol_main.jl, not just a slower run.
#

@testitem "example: mixture (water + urea)" setup=[TestFiles] begin
    using PDBToolsBase: read_pdb

    input_file = TestFiles.test_file("input_files/mixture.inp")
    sys = Packmol.read_packmol_input(input_file)
    @test sys.nmols == 1400

    sys.output_file = tempname() * ".pdb"
    result = Packmol.packmol(sys; nloop=20)
    @test result === sys
    @test sys.status == :packing_ready

    @test isfile(sys.output_file)
    @test length(read_pdb(sys.output_file)) == length(sys.atoms)
    rm(sys.output_file; force=true)
end

@testitem "example: spherical (double-layered vesicle)" setup=[TestFiles] begin
    using PDBToolsBase: read_pdb

    input_file = TestFiles.test_file("input_files/spherical.inp")
    sys = Packmol.read_packmol_input(input_file)
    @test sys.nmols == 18234

    sys.output_file = tempname() * ".pdb"
    result = Packmol.packmol(sys; nloop=20)
    @test result === sys
    @test sys.status == :packing_ready

    @test isfile(sys.output_file)
    @test length(read_pdb(sys.output_file)) == length(sys.atoms)
    rm(sys.output_file; force=true)
end

@testitem "example: solvprotein (protein + water + ions)" setup=[TestFiles] begin
    using PDBToolsBase: read_pdb

    input_file = TestFiles.test_file("input_files/solvprotein.inp")
    sys = Packmol.read_packmol_input(input_file)
    @test sys.nmols == 16551

    sys.output_file = tempname() * ".pdb"
    result = Packmol.packmol(sys; nloop=20)
    @test result === sys
    @test sys.status == :packing_ready

    @test isfile(sys.output_file)
    @test length(read_pdb(sys.output_file)) == length(sys.atoms)
    rm(sys.output_file; force=true)
end

@testitem "example: bilayer (lipid double layer)" setup=[TestFiles] begin
    using PDBToolsBase: read_pdb

    input_file = TestFiles.test_file("input_files/bilayer.inp")
    sys = Packmol.read_packmol_input(input_file)
    @test sys.nmols == 1100

    sys.output_file = tempname() * ".pdb"
    result = Packmol.packmol(sys; nloop=50)
    @test result === sys
    @test sys.status == :packing_ready

    @test isfile(sys.output_file)
    @test length(read_pdb(sys.output_file)) == length(sys.atoms)
    rm(sys.output_file; force=true)
end

@testitem "Fortran example: interface (xyz files)" setup=[TestFiles] begin
    using PDBToolsBase: read_xyz
    # Fortran Packmol's examples/interface.inp: xyz input and output, and a fixed
    # structure without `number`
    input_file = TestFiles.test_file("fortran_examples/interface.inp")
    sys = redirect_stdout(devnull) do
        Packmol.packmol(input_file; iprint=1000)
    end
    @test sys.status == :packing_ready
    output = read_xyz(joinpath(dirname(input_file), "interface.xyz"))
    @test length(output) == 1019 * 3 + 199 * 5 + 35
end
