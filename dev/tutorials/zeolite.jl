#
# Images of the zeolite tutorial (docs/src/tutorials/zeolite.md): packs the input of
# docs/src/tutorials/files/zeolite with the executable (juliac/build.jl), and renders it.
#
#     julia dev/tutorials/zeolite.jl
#
include("render.jl")
const exe = joinpath(@__DIR__, "..", "..", "juliac", "build", "bin", "packmol")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const work = mktempdir()
cp(joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files", "zeolite"), work; force=true)
run(pipeline(`$exe -i $(joinpath(work, "zeolite.inp"))`; stdout=devnull))
# Down the c axis of the hexagonal cell, slightly tilted
render(joinpath(work, "zeolite.pdb"), joinpath(images, "zeolite.webp");
    reps = [
        (selection = "resname CHA", style = "Licorice 0.25 12 12", color = "Element"),
        (selection = "resname CH4", style = "VDW 0.8 20", color = "ColorID 10"),
        (selection = "resname CO2", style = "VDW 0.8 20", color = "ColorID 7"),
    ],
    rotate = [(:x, -20), (:y, 10)], cell = :origin,
)
