module CLI

import ..Packmol

const usage = """
    Usage: packmol input.inp
           packmol -i input.inp

    With no input file, a file dialog is opened to choose it.
    """

function (@main)(ARGS)
    # Print to the `Core` streams (see `src/io/terminal.jl`)
    Packmol._use_core_streams[] = true
    if isempty(ARGS)
        input_file = Packmol.pick_file()
    elseif length(ARGS) == 1 && !startswith(ARGS[1], "-")
        input_file = ARGS[1]
    elseif length(ARGS) == 2 && ARGS[1] == "-i"
        input_file = ARGS[2]
    else
        Packmol._print(usage)
        return 1
    end
    # Dimension and number type given as types, so that the parser is type-stable
    system = Packmol.packmol(Packmol.read_packmol_input(input_file, Val(3), Float64))
    return system.status == :failed_packing ? 1 : 0
end

end
