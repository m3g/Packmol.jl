#
# Images of the lipid tutorials (POPC bilayer and vesicle): the packed system, and the
# initial approximation written with `check`. The inputs are those of
# docs/src/tutorials/files/{bilayer,vesicle}. Run from the root of the repository:
#
#     julia --project=docs -t 8 dev/tutorials/lipids.jl
#
using Packmol
include("render.jl")

const files = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const work = mktempdir()

# Runs the example (and its `check` version) and returns the output files
function pack(name, output)
    for f in readdir(joinpath(files, name))
        cp(joinpath(files, name, f), joinpath(work, f); force=true)
    end
    input = read(joinpath(work, "$name.inp"), String)
    write(joinpath(work, "check_$name.inp"), replace(input, "output $output" => "output check_$output\ncheck"))
    packmol(joinpath(work, "check_$name.inp"))
    packmol(joinpath(work, "$name.inp"))
    return joinpath(work, output), joinpath(work, "check_$output")
end

function both(name, output, image; kargs...)
    final, check = pack(name, output)
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

# A slab through the center of the vesicle, to show the layers
both("vesicle", "vesicle.pdb", "spherical";
    reps = [
        (selection = "water and same residue as (abs(y) < 3)", style = "Licorice 0.2 12 12"),
        (selection = "not water and same residue as (abs(y) < 3)", style = "Licorice 0.3 12 12"),
    ],
    rotate = [(:x, -90)],
)
