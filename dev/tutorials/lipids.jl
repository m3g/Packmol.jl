#
# Images of the lipid tutorials (POPC bilayer and nanoliposomes): the packed system, and the
# initial placement, written with the `check_initial_placement` line of the input
# uncommented. The inputs are those of
# docs/src/tutorials/files/{bilayer,nanoliposome}. Run from the root of the repository:
#
#     julia --project=docs -t 14 dev/tutorials/lipids.jl
#
using Packmol
include("render.jl")

const files = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const work = mktempdir()

# Runs the example `input` of the files of `name` (and its check version, unless
# `check = false`) and returns the output files
function pack(name, output; input = name, check = true)
    for f in readdir(joinpath(files, name))
        cp(joinpath(files, name, f), joinpath(work, f); force=true)
    end
    if check
        text = read(joinpath(work, "$input.inp"), String)
        write(joinpath(work, "check_$input.inp"), replace(text,
            "# check_initial_placement" => "check_initial_placement", "output $output" => "output check_$output"))
        packmol(joinpath(work, "check_$input.inp"))
    end
    packmol(joinpath(work, "$input.inp"))
    return joinpath(work, output), joinpath(work, "check_$output")
end

function both(name, output, image; input = name, kargs...)
    final, check = pack(name, output; input)
    render(final, joinpath(images, "$image.webp"); kargs...)
    render(check, joinpath(images, "$(image)_check.webp"); kargs...)
end

# The bilayer seen from the side, with the edges of the periodic cell
both("bilayer", "bilayer.pdb", "bilayer";
    reps = [
        (selection = "water", style = "Licorice 0.2 12 12"),
        (selection = "not water", style = "Licorice 0.3 12 12"),
    ],
    rotate = [(:x, -90), (:y, 20), (:x, 10)],
    cell = :origin,
)

# A slab through the center of the nanoliposomes, to show the layers (the larger one
# without the image of its initial approximation, which the tutorial does not show)
const slab = (
    reps = [
        (selection = "water and same residue as (abs(y) < 3)", style = "Licorice 0.2 12 12"),
        (selection = "not water and same residue as (abs(y) < 3)", style = "Licorice 0.3 12 12"),
    ],
    rotate = [(:x, -90)],
)
both("nanoliposome", "small.pdb", "nanoliposome_small"; input = "small", slab...)
final, _ = pack("nanoliposome", "large.pdb"; input = "large", check = false)
render(final, joinpath(images, "nanoliposome_large.webp"); slab...)
