#
# Regression tests based on the classic Fortran Packmol example inputs
# (full molecule counts, not the scaled-down variants used elsewhere in
# test/input_files). Each one must actually converge, and in a number of
# packing loops not much greater than what it currently needs — a large
# increase would signal a regression in the stall-detection/movebad!
# machinery in packmol_main.jl, not just a slower run.
#

@testitem "example: mixture (water + urea)" begin
    using PDBTools: read_pdb

    input_file = joinpath(Packmol.src_dir, "..", "test", "input_files", "mixture.inp")
    sys = Packmol.read_packmol_input(input_file)
    @test sys.nmols == 1400

    sys.output_file = tempname() * ".pdb"
    converged = Packmol.packmol(sys; nloop=20)
    @test converged

    @test isfile(sys.output_file)
    @test length(read_pdb(sys.output_file)) == length(sys.atoms)
    rm(sys.output_file; force=true)
end

@testitem "example: spherical (double-layered vesicle)" begin
    using PDBTools: read_pdb

    input_file = joinpath(Packmol.src_dir, "..", "test", "input_files", "spherical.inp")
    sys = Packmol.read_packmol_input(input_file)
    @test sys.nmols == 18234

    sys.output_file = tempname() * ".pdb"
    converged = Packmol.packmol(sys; nloop=20)
    @test converged

    @test isfile(sys.output_file)
    @test length(read_pdb(sys.output_file)) == length(sys.atoms)
    rm(sys.output_file; force=true)
end

@testitem "example: solvprotein (protein + water + ions)" begin
    using PDBTools: read_pdb

    input_file = joinpath(Packmol.src_dir, "..", "test", "input_files", "solvprotein.inp")
    sys = Packmol.read_packmol_input(input_file)
    @test sys.nmols == 16551

    sys.output_file = tempname() * ".pdb"
    converged = Packmol.packmol(sys; nloop=20)
    @test converged

    @test isfile(sys.output_file)
    @test length(read_pdb(sys.output_file)) == length(sys.atoms)
    rm(sys.output_file; force=true)
end

@testitem "example: bilayer (lipid double layer)" begin
    using PDBTools: read_pdb

    input_file = joinpath(Packmol.src_dir, "..", "test", "input_files", "bilayer.inp")
    sys = Packmol.read_packmol_input(input_file)
    @test sys.nmols == 1100

    sys.output_file = tempname() * ".pdb"
    converged = Packmol.packmol(sys; nloop=50)
    @test converged

    @test isfile(sys.output_file)
    @test length(read_pdb(sys.output_file)) == length(sys.atoms)
    rm(sys.output_file; force=true)
end
