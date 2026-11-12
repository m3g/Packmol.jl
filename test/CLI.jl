@testitem "CLI" setup=[TestFiles] begin
    import Packmol.CLI
    water = TestFiles.test_file("structure_files/water.pdb")
    dir = mktempdir()
    input = joinpath(dir, "box.inp")
    output = joinpath(dir, "other.pdb")
    write(input, """
    tolerance 2.0
    output box.pdb
    structure $water
      number 10
      inside box 0. 0. 0. 20. 20. 20.
    end structure
    """)
    failing_input = joinpath(dir, "failed.inp")
    write(failing_input, """
    tolerance 2.0
    output failed.pdb
    nloop 1
    maxit 2
    structure $water
      number 2
      inside box 0. 0. 0. 1.1 1.1 1.1
    end structure
    """)
    redirect_stdout(devnull) do
        @test CLI.main(["-h"]) == 0
        @test CLI.main(["--version"]) == 0
        @test CLI.main(["-x"]) == CLI.exit_code_command_line_error
        @test CLI.main(["a.inp", "b.inp"]) == CLI.exit_code_command_line_error
        @test CLI.main([joinpath(dir, "missing.inp")]) == CLI.exit_code_open_file
        @test CLI.main([input]) == CLI.exit_code_success
        @test isfile(joinpath(dir, "box.pdb"))
        # -o overrides the output of the input file, in any argument order
        @test CLI.main(["-o", output, "-i", input]) == CLI.exit_code_success
        @test isfile(output)
        # -s overrides the seed of the input file: the same seed gives the same
        # configuration, another seed a different one
        outputs = [joinpath(dir, "seed_$i.pdb") for i in 1:3]
        for (out, seed) in zip(outputs, ["7", "7", "8"])
            @test CLI.main(["-i", input, "-o", out, "--seed", seed]) == CLI.exit_code_success
        end
        @test read(outputs[1], String) == read(outputs[2], String)
        @test read(outputs[1], String) != read(outputs[3], String)
        @test CLI.main(["-i", input, "-s", "abc"]) == CLI.exit_code_command_line_error
        @test CLI.main([failing_input]) == CLI.exit_code_failed_to_converge
        @test isfile(joinpath(dir, "failed-FORCED.pdb"))
        write(input, "tolerance abc\n")
        @test CLI.main([input]) == CLI.exit_code_input_error
    end
    # `CLI.main` switches the printing to the `Core` streams
    Packmol._use_core_streams[] = false
end
