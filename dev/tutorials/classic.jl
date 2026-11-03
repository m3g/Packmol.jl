#
# Images of the classic tutorials (the examples of the Packmol site, with the files of
# docs/src/tutorials/files): the packed system, and the initial placement, written with the
# `check_initial_placement` line of the input uncommented (the bilayer and the nanoliposomes
# are in lipids.jl). Run from the root of the repository, with the packmol executable built
# (juliac/build.jl):
#
#     julia --project=docs dev/tutorials/classic.jl [example ...]
#
# (all the examples, or those given: mixture, interface, solvprotein)
#
include("render.jl")

const exe = joinpath(@__DIR__, "..", "..", "juliac", "build", "bin", "packmol")
const files = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const work = mktempdir()

# The input with the `check_initial_placement` line uncommented, writing to `check_<output>`
check_input(input, output) = replace(input,
    "# check_initial_placement" => "check_initial_placement",
    "output $output" => "output check_$output",
)

# Runs the example (and its check version) and returns the output files
function pack(name, output)
    for f in readdir(joinpath(files, name))
        cp(joinpath(files, name, f), joinpath(work, f); force=true)
    end
    input = read(joinpath(work, "$name.inp"), String)
    write(joinpath(work, "check_$name.inp"), check_input(input, output))
    run(pipeline(`$exe -i $(joinpath(work, "check_$name.inp"))`; stdout=devnull))
    run(pipeline(`$exe -i $(joinpath(work, "$name.inp"))`; stdout=devnull))
    return joinpath(work, output), joinpath(work, "check_$output")
end

function both(name, output; check_image=true, kargs...)
    isempty(ARGS) || name in ARGS || return nothing
    final, check = pack(name, output)
    render(final, joinpath(images, "$name.webp"); kargs...)
    check_image && render(check, joinpath(images, "$(name)_check.webp"); kargs...)
end

# A simple system: no `check` image
both("mixture", "mixture.pdb"; check_image=false,
    reps = [
        (selection = "water", style = "Licorice 0.2 12 12"),
        (selection = "not water", style = "Licorice 0.35 12 12"),
    ],
    rotate = [(:y, 30), (:x, 20)],
)

# Water (left), chloroform (right, as spheres: the XYZ files have no residues, so bonds
# would be guessed between molecules) and the hormone at the interface, on top of the
# solvent below its plane (z < 20)
both("interface", "interface.xyz";
    reps = [
        (selection = "index < 3057 and z < 20", style = "Licorice 0.2 12 12"),
        (selection = "index >= 3057 and index < 4052 and z < 20", style = "VDW 0.4 12", color = "ColorID 7"),
        (selection = "index >= 4052", style = "VDW 1.0 20", material = "Glossy"),
    ],
    rotate = [(:x, -30)],
)

# The water behind the protein (y > 0, away from the viewer), sodium (orange) and chloride
# (green), in the rhombic dodecahedron (seen along y, a hexagon)
both("solvprotein", "solvprotein.pdb";
    reps = [
        (selection = "protein", style = "NewCartoon 0.3 20 4.1 0", color = "ColorID 0"),
        (selection = "water and same residue as (y > 0)", style = "Lines 1"),
        (selection = "resname SOD", style = "VDW 1.0 20", color = "ColorID 3"),
        (selection = "resname CLA", style = "VDW 1.0 20", color = "ColorID 7"),
    ],
    rotate = [(:x, -90)],
)
