#
# Images of the Zika virion tutorial (docs/src/tutorials/zika.md), from the outputs of
# the input of the packmol_zika set of MolSimTestFiles (with and without `check`), given
# as arguments: `julia dev/tutorials/zika.jl zika_check.pdb zika.pdb`.
#
include("render.jl")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
# A slab through the center of the virion (|y| < 10 Å), seen along y
slab = "same residue as (abs(y) < 10)"
reps = [
    (selection = "index < 489960 and abs(y) < 10", style = "VDW 1.6 12", color = "ColorID 0"),
    (selection = "resname POC POE POS and $slab", style = "VDW 1.4 12", color = "ColorID 3"),
    (selection = "resname WT4 WLS and $slab", style = "VDW 1.0 8", color = "ColorID 15"),
    (selection = "resname NaW ClW and $slab", style = "VDW 1.6 12", color = "ResName"),
]
for (file, name) in zip(ARGS, ("zika_check.webp", "zika.webp"))
    render(file, joinpath(images, name); reps, rotate = [(:x, -90)])
end
