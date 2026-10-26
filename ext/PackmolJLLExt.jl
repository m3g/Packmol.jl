module PackmolJLLExt

import Packmol
import Packmol_jll

function Packmol.run_packmol_legacy(input_file::String)
    input_file = expanduser(input_file)
    isfile(input_file) || error("Input file not found: $input_file")
    current_dir = pwd()
    println("Current directory: $current_dir")
    input_file_dir = dirname(abspath(input_file))
    println("Input file directory: $input_file_dir")
    output_file_name = normpath(joinpath(input_file_dir,
        open(input_file) do io
            for line in eachline(io)
                if occursin("output", line)
                    _, output_file = split(line)
                    return output_file
                end
            end
        end
    ))
    println("Output file: $output_file_name")
    # Run packmol
    try
        cd(input_file_dir)
        run(pipeline(`$(Packmol_jll.packmol())`, stdin=input_file))
        println("Wrote output to: ", output_file_name)
        cd(current_dir)
    catch
        cd(current_dir)
        @error "Error running packmol" _file=nothing _line=nothing
    end
    return nothing
end

Packmol.run_packmol_legacy() = Packmol.run_packmol_legacy(Packmol.pick_file())

end
