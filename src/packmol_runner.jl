import Packmol_jll
import NativeFileDialog

packmol_runner = Packmol_jll.packmol()

export run_packmol

function run_packmol(input_file::String)
    isfile(input_file) || error("Input file not found: $input_file")
    cd(dirname(input_file))
    # Run packmol  
    run(pipeline(`$packmol_runner`, stdin=input_file))
    println("Wrote output to: ", dirname(input_file))
    return nothing
end

function run_packmol()
    input_file = NativeFileDialog.pick_file() 
    run_packmol(input_file)
end

@testitem "run_packmol" begin
    run_packmol("../test/run_packmol/water_box.inp")
    @test isfile("../test/run_packmol/water_box.pdb")
    run_packmol("../test/run_packmol/ieee_signaling.inp")
    @test isfile("../test/run_packmol/ieee_signaling_box.pdb")
end


