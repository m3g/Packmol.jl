#
# Images of the protein in glycerol/water tutorial (docs/src/tutorials/protein_glycerol.md):
# packs the input of docs/src/tutorials/files/protein_glycerol with the executable
# (juliac/build.jl), and renders it.
#
#     julia dev/tutorials/protein_glycerol.jl
#
include("render.jl")
const exe = joinpath(@__DIR__, "..", "..", "juliac", "build", "bin", "packmol")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const work = mktempdir()
cp(joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files", "protein_glycerol"), work; force=true)
run(pipeline(`$exe -i $(joinpath(work, "protein_glycerol.inp"))`; stdout=devnull))
# The solvent behind the protein (y > 0, away from the viewer)
render(joinpath(work, "protein_glycerol.pdb"), joinpath(images, "protein_glycerol.webp");
    reps = [
        (selection = "protein", style = "NewCartoon 0.3 20 4.1 0", color = "ColorID 0"),
        (selection = "resname GLYC and same residue as (y > 0)", style = "Licorice 0.3 12 12", color = "ColorID 7"),
        (selection = "water and same residue as (y > 0)", style = "Lines 1"),
    ],
    rotate = [(:x, -90)], cell = :origin,
)
