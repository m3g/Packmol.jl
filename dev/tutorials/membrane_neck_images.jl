#
# Image of the membrane tube with a neck tutorial: packs the system with the script of
# docs/src/tutorials/files/membrane_neck and renders a slab through the axis of the tube.
# Run from the root of the repository:
#
#     julia --project=docs -t 14 dev/tutorials/membrane_neck_images.jl
#
include("render.jl")

const files = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files", "membrane_neck")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const work = mktempdir()

using Packmol, LazyArtifacts # Packmol.tutorial_file
cp(joinpath(files, "membrane_neck.jl"), joinpath(work, "membrane_neck.jl"); force = true)
for f in ("popc.pdb", "water.pdb")
    cp(Packmol.tutorial_file(f), joinpath(work, f); force = true)
end
cd(() -> include(joinpath(work, "membrane_neck.jl")), work)

render(joinpath(work, "membrane_neck.pdb"), joinpath(images, "membrane_neck.webp");
    reps = [
        (selection = "water and same residue as (abs(x) < 3)", style = "Licorice 0.2 12 12"),
        (selection = "not water and same residue as (abs(x) < 3)", style = "Licorice 0.3 12 12"),
    ],
    rotate = [(:y, 90)], size = (2000, 1400), width = 1400)
