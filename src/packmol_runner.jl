import Packmol_jll

packmol_runner = Packmol_jll.packmol()

export run_packmol

"""
    run_packmol()
    run_packmol(input_file::String)

Runs the packmol executable with the input file `input_file`. This will run
the classical `http://m3g.iqm.unicamp.br/packmol` program, which is a
pre-compiled binary. The input file is a text file with the same syntax as
the packmol input files.

If no input file is provided, a file explorer will be opened to choose the 
input file.

"""
function run_packmol end

@doc (@doc run_packmol)
function run_packmol(input_file::String)
    isfile(input_file) || error("Input file not found: $input_file")
    current_dir = pwd()
    println("Current directory: $current_dir")
    input_file_dir = dirname(abspath(input_file))
    println("Input file directory: $input_file_dir")
    ouput_file_name = normpath(joinpath(input_file_dir,
        open(input_file) do io
            for line in eachline(io)
                if occursin("output", line)
                    _, output_file = split(line)
                    return output_file
                end
            end
        end
    ))
    println("Output file: $ouput_file_name")
    # Run packmol  
    try 
        cd(input_file_dir)
        run(pipeline(`$packmol_runner`, stdin=input_file))
        println("Wrote output to: ", ouput_file_name)
        cd(current_dir)
    catch 
        cd(current_dir)
        @error "Error running packmol" _file=nothing _line=nothing
    end
    return nothing
end

@testitem "run_packmol" begin
    using Packmol
    test_dir = Packmol.src_dir*"/../test/run_packmol"
    cd(test_dir)
    run_packmol("water_box.inp")
    @test isfile("water_box.pdb")
    run_packmol("ieee_signaling.inp")
    @test isfile("ieee_signaling_box.pdb")
    cd(test_dir)
    rm("water_box.pdb")
    run_packmol("water_box.inp")
    @test isfile("water_box.pdb")
end

@static if haskey(ENV, "PACKMOL_GUI") && ENV["PACKMOL_GUI"] == "false"
    function run_packmol()
        throw(ArgumentError("""\n
            Environment variable PACKMOL_GUI is set to false. Set it to true to use the file dialog.
            
        """))
    end
else
    import NativeFileDialog
    @doc (@doc run_packmol)
    function run_packmol()
        input_file = NativeFileDialog.pick_file() 
        run_packmol(input_file)
    end
end

