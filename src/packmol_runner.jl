import Packmol_jll

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


