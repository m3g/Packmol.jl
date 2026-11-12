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

    @testset "output redirected; screen_log" begin
        # With stdout redirected (`packmol input.inp > log`), nothing is printed to stderr,
        # and the log has the report of the run, with the progress bars (no carriage
        # returns or escape sequences)
        log = joinpath(dir, "redirected.log")
        err = IOBuffer()
        cmd = Cmd(`$exe -i $(joinpath(inputs, "water_box.inp"))`; ignorestatus=true, dir=inputs)
        p = run(pipeline(cmd; stdout=log, stderr=err))
        @test p.exitcode == 0
        @test isempty(take!(err))
        report = read(log, String)
        @test occursin("Success!", report)
        @test occursin("|*****", report)
        @test !occursin('\r', report) && !occursin('\e', report)
        # `screen_log no`: nothing printed, also not the warnings
        quiet = joinpath(inputs, "water_box_quiet.inp")
        write(quiet, read(joinpath(inputs, "water_box.inp"), String) * "\nfscale 1.0\nscreen_log no\n")
        out = IOBuffer()
        err = IOBuffer()
        p = run(pipeline(Cmd(`$exe -i $quiet`; ignorestatus=true, dir=inputs); stdout=out, stderr=err))
        @test p.exitcode == 0
        @test isempty(take!(out))
        @test isempty(take!(err))
        # Errors are printed to stderr
        bad = joinpath(dir, "missing_structure.inp")
        write(bad, "tolerance 2.0\noutput x.pdb\nstructure missing.pdb\n  number 1\n  inside cube 0. 0. 0. 10.\nend structure\n")
        out = IOBuffer()
        err = IOBuffer()
        p = run(pipeline(Cmd(`$exe -i $bad`; ignorestatus=true, dir=dir); stdout=out, stderr=err))
        @test p.exitcode == 172
        @test occursin("ERROR: opening file", String(take!(err)))
        @test !occursin("ERROR", String(take!(out)))
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

    @testset "mmCIF output" begin
        input_file = joinpath(inputs, "water_box_pbc.inp")
        output_file = joinpath(dir, "water_box_pbc.cif")
        code, out = run_packmol("-i", input_file, "-o", output_file)
        @test code == 0
        @test isfile(output_file)
        atoms = Packmol.PDBToolsBase.read_mmcif(output_file)
        @test length(atoms) == length(read_packmol_input(input_file).atoms)
        @test Packmol.PDBToolsBase.read_unitcell(output_file) ≈ read_packmol_input(input_file).unitcell rtol = 1e-3
    end
    @testset "volume estimates" begin
        # The executable gives the estimates of the Julia functions (the same random points)
        for name in ("water_box", "solvprotein_pbc")
            input_file = joinpath(inputs, "$name.inp")
            code, out = run_packmol("--guess-volume", "--density", "1.1", input_file)
            @test code == 0
            g = redirect_stdout(devnull) do
                Packmol.guess_volume(input_file; density=1.1)
            end
            @test occursin("Volume:       $(Packmol._fmt(g.volume, 1)) Å³", out)
            code, out = run_packmol("--free-volume", "--npoints", "20000", input_file)
            @test code == 0
            r = redirect_stdout(devnull) do
                Packmol.free_volume(input_file; npoints=20000)
            end
            @test all(x -> occursin(Packmol._fmt(x.volume, 1), out), r)
        end
        code, _ = run_packmol("--free-volume", "--npoints", "0", joinpath(inputs, "water_box.inp"))
        @test code == 174
        code, _ = run_packmol("--guess-volume", "--density", "-1", joinpath(inputs, "water_box.inp"))
        @test code == 174
    end
    @testset "voxel_place" begin
        # The voxel placement of Bentopy, then the packing of the molecules that fit
        vdir = mktempdir()
        cp(joinpath(dir, "structure_files", "water.pdb"), joinpath(vdir, "water.pdb"))
        cp(joinpath(dir, "structure_files", "protein.pdb"), joinpath(vdir, "protein.pdb"))
        input_file = joinpath(vdir, "voxel.inp")
        write(input_file, """
            tolerance 2.0
            voxel_place 500
            pbc 120. 120. 120.
            output voxel.pdb
            structure protein.pdb
              number 1
              fixed 60. 60. 60. 0. 0. 0.
            end structure
            structure water.pdb
              number 4000
            end structure
            """)
        code, out = run_packmol("-i", input_file)
        @test code == 0
        @test occursin("Placing the molecules on a voxel grid", out)
        @test occursin("Placed 4000 of 4000 free molecules", out)
        dmin, sys = minimum_intermolecular_distance(input_file)
        @test dmin >= sys.tolerance - sys.tolerance_precision - pdb_rounding
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

    @testset "surface constraints" begin
        # Water within a shell around a sphere given by points with normals, the oxygen
        # atoms not within 4 Å of the sphere
        sdir = mktempdir()
        cp(joinpath(dir, "structure_files", "water.pdb"), joinpath(sdir, "water.pdb"))
        open(joinpath(sdir, "sphere.dat"), "w") do io
            n = 1500
            for i in 1:n
                θ = acos(1 - 2 * (i - 0.5) / n); φ = π * (1 + sqrt(5)) * i
                u = (sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))
                println(io, join(15 .* u, " "), " ", join(u, " "))
            end
        end
        input_file = joinpath(sdir, "shell.inp")
        write(input_file, """
            tolerance 2.0
            output shell.pdb
            structure water.pdb
              number 120
              within surface sphere.dat 2.0 7.0
              atoms 1
                not within surface sphere.dat 4.0
              end atoms
            end structure
            """)
        code, out = run_packmol("-i", input_file; cwd=sdir)
        @test code == 0
        code == 0 || println(out)
        r = [sqrt(sum(abs2, position(at))) for at in read_pdb(joinpath(sdir, "shell.pdb"))]
        @test 19.5 < minimum(r) && maximum(r) < 24.5
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
        # Failed packing: the forced solution is written to <name>-FORCED.<ext> (Fortran:
        # <output>_FORCED, test_failed.sh).
        # The box of water_box_failed.inp (3 x 3 x 1 Å) holds the two molecules (they fail
        # in Fortran with its `maxit 2`), so it is replaced by one that cannot: its
        # diagonal (1.91 Å) is shorter than the tolerance
        failed = joinpath(inputs, "water_box_failed.inp")
        input_file = joinpath(inputs, "water_box_impossible.inp")
        write(input_file, replace(read(failed, String), "inside box 0. 0. 0. 3.0 3.0 1.0" => "inside box 0. 0. 0. 1.1 1.1 1.1"))
        code, out = run_packmol("-i", input_file)
        @test code == 173
        @test occursin("FORCED", out)
        @test isfile(Packmol._forced_file(read_packmol_input(input_file).output_file))
        # Constraints that the initial approximation cannot satisfy: an input error,
        # unless continue_failed_initial_packing is set
        bad = joinpath(dir, "bad_constraints.inp")
        write(bad, """
            tolerance 2.0
            output $(joinpath(dir, "bad_constraints.pdb"))
            nloop 5
            structure $(joinpath(inputs, "water.pdb"))
              number 10
              inside box 0. 0. 0. 10. 10. 10.
              outside box -1. -1. -1. 11. 11. 11.
            end structure
            """)
        code, out = run_packmol("-i", bad)
        @test code == 171
        @test occursin("continue_failed_initial_packing", out)
        write(bad, replace(read(bad, String), "nloop 5" => "nloop 5\ncontinue_failed_initial_packing"))
        code, _ = run_packmol("-i", bad)
        @test code == 173
        # Input error
        code, out = run_packmol("-i", joinpath(inputs, "protein_outside_pbc_error.inp"))
        @test code == 171
        @test occursin("ERROR", out)
        # Input file not found
        code, _ = run_packmol("-i", joinpath(inputs, "does_not_exist.inp"))
        @test code in (171, 172)
    end
end
