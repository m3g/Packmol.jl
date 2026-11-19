#
# Images of the condensate tutorial (docs/src/tutorials/condensate.md), from the outputs
# of its input files (docs/src/tutorials/files/condensate), packed in that directory:
#
#     packmol homogeneous.inp
#     packmol droplet_chains.inp
#     packmol droplet_water.inp
#
# and the initial placement of droplet_chains.inp, written to droplet_check.pdb:
#
#     julia --project -e 'using Packmol; packmol("droplet_chains.inp"; check_initial_placement = true, output = "droplet_check.pdb")'
#
# then, from the root of the repository:
#
#     julia --project=dev/large_tests dev/tutorials/condensate_images.jl
#
# The chains are colored by chain; one water bead (WN1) in 20 is shown in the whole
# views, and all of them in the slabs (|z| < 10 Å).
#
include("render.jl")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const dir = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files", "condensate")

const chains = (selection = "resname sV sP sG", style = "VDW 1.6 12", color = "Chain")
const some_water = (selection = "resname WT4 and name WN1 and resid % 20 == 0", style = "VDW 1.6 12",
    color = "ColorID 10", material = "Transparent")
slab(rep) = merge(rep, (selection = "($(rep.selection)) and abs(z) < 10",))
const all_water = (selection = "resname WT4", style = "VDW 1.6 12", color = "ColorID 10", material = "Transparent")

render(joinpath(dir, "droplet_check.pdb"), joinpath(images, "condensate_check.webp");
    reps = [chains], rotate = [(:z, 30), (:x, -70)], cell = :origin)
render(joinpath(dir, "droplet.pdb"), joinpath(images, "condensate_droplet.webp");
    reps = [chains, some_water], rotate = [(:z, 30), (:x, -70)], cell = :origin)
render(joinpath(dir, "droplet.pdb"), joinpath(images, "condensate_droplet_slab.webp");
    reps = [slab(chains), slab(all_water)], cell = :origin)
render(joinpath(dir, "homogeneous.pdb"), joinpath(images, "condensate_homogeneous.webp");
    reps = [chains, some_water], rotate = [(:z, 30), (:x, -70)], cell = :origin)
render(joinpath(dir, "homogeneous.pdb"), joinpath(images, "condensate_homogeneous_slab.webp");
    reps = [slab(chains), slab(all_water)], cell = :origin)
