#
# Images of the tutorial of the polymer-capped gold nanoparticle
# (docs/src/tutorials/gold_polymer.md), from the outputs of its input file
# (docs/src/tutorials/files/gold_polymer), packed in that directory:
#
#     packmol gold_polymer.inp
#
# and its initial placement, written to gold_polymer_check.pdb:
#
#     julia --project -e 'using Packmol; packmol("gold_polymer.inp"; check_initial_placement = true, output = "gold_polymer_check.pdb")'
#
# then, from the root of the repository:
#
#     julia --project=dev/large_tests dev/tutorials/gold_polymer_images.jl
#
# The gold core in yellow, the sulfur atoms in orange, the chains colored by chain; in the
# whole views, one in five toluene molecules of the half of the cell behind the particle,
# and all of them in the slab (|z| < 5 Å).
#
include("render.jl")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const dir = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files", "gold_polymer")

const gold = (selection = "resname AU", style = "VDW 1.0 12", color = "ColorID 4")
const sulfur = (selection = "resname ETS and element S", style = "VDW 1.0 12", color = "ColorID 3")
const chains = (selection = "resname ETS STY SBU", style = "Licorice 0.3 12 12", color = "Chain")
const some_toluene = (selection = "resname TOL and resid % 5 == 0 and y > 0", style = "Licorice 0.3 12 12",
    color = "ColorID 7", material = "Transparent")
slab(rep) = merge(rep, (selection = "($(rep.selection)) and abs(z) < 5",))
const all_toluene = (selection = "resname TOL", style = "Licorice 0.3 12 12", color = "ColorID 7")

render(joinpath(dir, "gold_polymer_check.pdb"), joinpath(images, "gold_polymer_check.webp");
    reps = [gold, sulfur, chains], rotate = [(:z, 30), (:x, -70)], cell = :origin)
render(joinpath(dir, "gold_polymer.pdb"), joinpath(images, "gold_polymer.webp");
    reps = [gold, sulfur, chains, some_toluene], rotate = [(:z, 30), (:x, -70)], cell = :origin)
render(joinpath(dir, "gold_polymer.pdb"), joinpath(images, "gold_polymer_slab.webp");
    reps = [slab(gold), slab(sulfur), slab(chains), slab(all_toluene)], cell = :origin)
