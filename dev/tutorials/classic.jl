#
# Images of the classic tutorials (the examples of the Packmol site): the packed system,
# and the initial approximation written with `check`. Run from the root of the
# repository, with the packmol executable built (juliac/build.jl):
#
#     julia --project=docs dev/tutorials/classic.jl
#
using Packmol, LazyArtifacts
include("render.jl")

const exe = joinpath(@__DIR__, "..", "..", "juliac", "build", "bin", "packmol")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const work = mktempdir()
cp(Packmol.test_file("fortran_examples"), work; force=true)
chmod(work, 0o755; recursive=true)

# Runs the example (and its `check` version) and returns the output files
function pack(name, output)
    input = read(joinpath(work, "$name-comment.inp"), String)
    check = replace(input, "output $output" => "output check_$output\ncheck")
    write(joinpath(work, "check_$name.inp"), check)
    run(pipeline(`$exe -i $(joinpath(work, "check_$name.inp"))`; stdout=devnull))
    run(pipeline(`$exe -i $(joinpath(work, "$name-comment.inp"))`; stdout=devnull))
    return joinpath(work, output), joinpath(work, "check_$output")
end

function both(name, output; check_image=true, kargs...)
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

both("bilayer", "bilayer.pdb";
    reps = [
        (selection = "water", style = "Licorice 0.2 12 12"),
        (selection = "not water", style = "Licorice 0.3 12 12"),
    ],
    rotate = [(:x, -90), (:y, 20), (:x, 10)],
)

# A slab through the center of the vesicle, to show the layers
both("spherical", "spherical.pdb";
    reps = [
        (selection = "water and same residue as (abs(y) < 3)", style = "Licorice 0.2 12 12"),
        (selection = "not water and same residue as (abs(y) < 3)", style = "Licorice 0.3 12 12"),
    ],
    rotate = [(:x, -90)],
)

# The water behind the protein (y > 0, away from the viewer)
both("solvprotein", "solvprotein.pdb";
    reps = [
        (selection = "protein", style = "NewCartoon 0.3 20 4.1 0", color = "ColorID 0"),
        (selection = "water and same residue as (y > 0)", style = "Lines 1"),
        (selection = "resname SOD CLA", style = "VDW 1.0 20", color = "ResName"),
    ],
    rotate = [(:x, -90)],
)
