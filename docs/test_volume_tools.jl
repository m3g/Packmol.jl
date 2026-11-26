# The Julia estimates of free_volume on the inputs of the tutorials, for the tests of the
# JavaScript version (docs/test_volume_tools.js):
#   julia --project docs/test_volume_tools.jl <results.json>
using Packmol
using LazyArtifacts: LazyArtifacts # Packmol.tutorial_file
const files = joinpath(@__DIR__, "src", "tutorials", "files")
# Each input is copied, with the structure files of the tutorials (the packmol_tutorials
# set of MolSimTestFiles), to a directory next to the results, where the JavaScript tests
# read them
const work = abspath(ARGS[1]) * "_inputs"
const inputs = [
    "mixture/mixture.inp", "solvprotein/solvprotein.inp", "bilayer/bilayer.inp",
    "nanoliposome/small.inp", "water_drop/drop.inp", "zeolite/zeolite.inp",
    "interface/interface.inp", "nanotube/nanotube.inp",
]
const npoints = 200_000
json(x::Real) = isnan(x) ? "\"NaN\"" : string(x)
open(ARGS[1], "w") do io
    println(io, "[")
    for (k, f) in enumerate(inputs)
        dir = joinpath(work, dirname(f))
        if !isdir(dir)
            mkpath(dir)
            for g in readdir(joinpath(files, dirname(f)))
                cp(joinpath(files, dirname(f), g), joinpath(dir, g))
            end
            for g in readdir(Packmol.tutorial_file(""))
                cp(Packmol.tutorial_file(g), joinpath(dir, g); force = true)
            end
        end
        input = joinpath(work, f)
        r = redirect_stdout(devnull) do
            free_volume(input; npoints)
        end
        items = join(["{\"volume\": $(json(x.volume)), \"error\": $(json(x.error)), \"density\": $(json(x.density))}" for x in r], ", ")
        # (the path escaped: the backslashes of Windows)
        println(io, "  {\"input\": \"$(escape_string(input))\", \"npoints\": $npoints, \"results\": [$items]}", k < length(inputs) ? "," : "")
    end
    println(io, "]")
end
