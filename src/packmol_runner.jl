import Packmol_jll
import NativeFileDialog

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
    cd(dirname(abspath(input_file)))
    # Run packmol  
    run(pipeline(`$packmol_runner`, stdin=input_file))
    println("Wrote output to: ", dirname(abspath(input_file)))
    return nothing
end

@doc (@doc run_packmol)
function run_packmol()
    input_file = NativeFileDialog.pick_file() 
    run_packmol(input_file)
end

@testitem "run_packmol" begin
    using Packmol
    test_dir = Packmol.src_dir*"/../test/run_packmol"
    run_packmol("$test_dir/water_box.inp")
    @test isfile("$test_dir/water_box.pdb")
    run_packmol("$test_dir/ieee_signaling.inp")
    @test isfile("$test_dir/ieee_signaling_box.pdb")
    cd(test_dir)
    rm("water_box.pdb")
    run_packmol("water_box.inp")
    @test isfile("water_box.pdb")
end


