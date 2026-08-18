module CLI

import Packmol_jll

function (@main)(ARGS)
    proc = run(`$(Packmol_jll.packmol()) $ARGS`)
    proc.exitcode
end

end