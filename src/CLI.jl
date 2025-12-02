module CLI

using ..Packmol_jll: packmol

function (@main)(ARGS)
    proc = run(`$(packmol()) $ARGS`)
    proc.exitcode
end

end