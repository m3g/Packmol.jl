# The Julia estimates of free_volume on the inputs of the tutorials, for the tests of the
# JavaScript version (docs/test_volume_tools.js):
#   julia --project docs/test_volume_tools.jl <results.json>
using Packmol
const files = joinpath(@__DIR__, "src", "tutorials", "files")
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
        input = joinpath(files, f)
        r = redirect_stdout(devnull) do
            free_volume(input; npoints)
        end
        items = join(["{\"volume\": $(json(x.volume)), \"error\": $(json(x.error)), \"density\": $(json(x.density))}" for x in r], ", ")
        println(io, "  {\"input\": \"$input\", \"npoints\": $npoints, \"results\": [$items]}", k < length(inputs) ? "," : "")
    end
    println(io, "]")
end
