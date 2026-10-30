#
# Tests of the compiled `packmol` executable (the drop-in acceptance test, after Fortran
# Packmol's testing/ scripts): every input of the test file set is packed by the
# executable, and the output is checked (atoms of different molecules at least
# `tolerance - tolerance_precision` apart), as are the command line interface and the
# exit codes.
#
# Usage, from the root of the repository (after juliac/build.jl):
#
#     julia --project=juliac juliac/test_executable.jl [path/to/packmol]
#
# The executable defaults to `juliac/build/bin/packmol`.
#
using Test
using LazyArtifacts: LazyArtifacts
using Packmol: Packmol, read_packmol_input
using PDBToolsBase: read_pdb, read_xyz, position
using StaticArrays: SVector
import CellListMap

const exe = abspath(isempty(ARGS) ? joinpath(@__DIR__, "build", "bin", Sys.iswindows() ? "packmol.exe" : "packmol") : ARGS[1])
isfile(exe) || error("executable not found: $exe")

# Writable copy of the test files: the executable writes the outputs next to the inputs
const dir = mktempdir()
cp(Packmol.test_file(""), dir; force=true)
chmod(dir, 0o755; recursive=true)
const inputs = joinpath(dir, "input_files")

# Maximum change of a distance by the rounding of the coordinates in the output files
# (0.001 Å in PDB files): the convergence criterion holds for the coordinates before rounding
const pdb_rounding = 2 * sqrt(3) * 0.0005

# Runs the executable, and returns its exit code and output (stdout and stderr)
function run_packmol(args...; stdin_file=nothing, cwd=inputs)
    cmd = Cmd(`$exe $args`; ignorestatus=true, dir=cwd)
    isnothing(stdin_file) || (cmd = pipeline(cmd; stdin=stdin_file))
    out = IOBuffer()
    p = run(pipeline(cmd; stdout=out, stderr=out))
    return p.exitcode, String(take!(out))
end

# Minimum distance between atoms of different molecules in the output of the packing
# described by the input file (the input file gives the molecules and the unit cell)
function minimum_intermolecular_distance(input_file, output_file=nothing)
    sys = read_packmol_input(input_file)
    output_file = something(output_file, sys.output_file)
    atoms = endswith(output_file, ".xyz") ? read_xyz(output_file) : read_pdb(output_file)
    @test length(atoms) == length(sys.atoms)
    x = [SVector{3,Float64}(position(at)) for at in atoms]
    unitcell = isnothing(sys.unitcell) ? nothing : sys.unitcell
    pairs = CellListMap.neighborlist(; positions=x, cutoff=1.5 * sys.tolerance, unitcell)
    dmin = Inf
    for (i, j, d) in pairs
        sys.atoms[i].molecule_index == sys.atoms[j].molecule_index && continue
        dmin = min(dmin, d)
    end
    return dmin, sys
end

@testset "packmol executable" begin
    @testset "version and help" begin
        code, out = run_packmol("--version")
        @test code == 0
        @test startswith(out, "packmol ")
        code, out = run_packmol("-h")
        @test code == 0
        @test occursin("-i", out)
        code, _ = run_packmol("--not-an-option")
        @test code == 174
        code, _ = run_packmol("-t", "0", "-i", joinpath(inputs, "water_box.inp"))
        @test code == 174
    end

    @testset "threads" begin
        input_file = joinpath(inputs, "water_box.inp")
        for (args, env, n) in [(("-t", "1"), nothing, 1), (("-t", "2"), "1", 2), ((), "2", 2)]
            cmd = Cmd(`$exe $args -i $input_file`; ignorestatus=true)
            cmd = addenv(cmd, "JULIA_NUM_THREADS" => something(env, ""))
            out = read(cmd, String)
            @test occursin(Regex("Number of threads:\\s+$n\\b"), out)
        end
    end

    @testset "$name" for name in [
        "water_box", "water_box_small", "water_box_conect", "ieee_signaling", "mixture",
        "mixture_pbc", "toy_with_fixed", "water_box_pbc", "water_box_pbc2",
        "water_box_pbc_negative_coordinates", "water_box_pbc_slab", "water_box_pbc_outside_box",
        "bilayer", "bilayer_pbc", "solvprotein", "solvprotein_pbc", "spherical", "spherical_pbc",
        "only_one_fixed", "benzene2",
    ]
        input_file = joinpath(inputs, "$name.inp")
        elapsed = @elapsed code, out = run_packmol("-i", input_file)
        @test code == 0
        code == 0 || println(out)
        dmin, sys = minimum_intermolecular_distance(input_file)
        @test dmin >= sys.tolerance - sys.tolerance_precision - pdb_rounding
        println("  $name: $(round(elapsed; digits=2)) s")
    end

    @testset "input from stdin, and -o" begin
        input_file = joinpath(inputs, "water_box.inp")
        code, _ = run_packmol(; stdin_file=input_file)
        @test code == 0
        dmin, sys = minimum_intermolecular_distance(input_file)
        @test dmin >= sys.tolerance - sys.tolerance_precision - pdb_rounding
        output_file = joinpath(dir, "water_box_o.pdb")
        code, _ = run_packmol("-i", input_file, "-o", output_file)
        @test code == 0
        @test isfile(output_file)
        dmin, sys = minimum_intermolecular_distance(input_file, output_file)
        @test dmin >= sys.tolerance - sys.tolerance_precision - pdb_rounding
        # The input file as the first positional argument
        rm(output_file)
        code, _ = run_packmol(input_file, "-o", output_file)
        @test code == 0
        @test isfile(output_file)
    end

    @testset "examples" begin
        examples = joinpath(dir, "fortran_examples")
        for name in ["interface"]
            input_file = joinpath(examples, "$name.inp")
            code, out = run_packmol("-i", input_file)
            @test code == 0
            code == 0 || println(out)
        end
    end

    @testset "exit codes" begin
        # Failed packing: the best solution is written to <output>_FORCED (test_failed.sh)
        input_file = joinpath(inputs, "water_box_failed.inp")
        code, out = run_packmol("-i", input_file)
        @test code == 173
        @test occursin("FORCED", out)
        @test isfile(read_packmol_input(input_file).output_file * "_FORCED")
        # Input error
        code, out = run_packmol("-i", joinpath(inputs, "protein_outside_pbc_error.inp"))
        @test code == 171
        @test occursin("ERROR", out)
        # Input file not found
        code, _ = run_packmol("-i", joinpath(inputs, "does_not_exist.inp"))
        @test code in (171, 172)
    end
end
