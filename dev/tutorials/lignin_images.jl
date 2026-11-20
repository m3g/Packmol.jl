#
# Images of the tutorial of the lignin particle (docs/src/tutorials/lignin.md), from the
# outputs of its input files (docs/src/tutorials/files/lignin), packed in that directory:
#
#     packmol particle.inp
#     packmol solvation.inp
#     packmol dispersed.inp
#
# and the initial placement of the particle, written to particle_check.pdb:
#
#     julia --project -e 'using Packmol; packmol("particle.inp"; check_initial_placement = true, output = "particle_check.pdb")'
#
# then, from the root of the repository:
#
#     julia --project=dev/large_tests dev/tutorials/lignin_images.jl
#
# The oligomers colored by molecule (chain); the solvent as transparent spheres, in the
# half of the cell beyond the cut x = 0 in the whole views, and in the slab (|z| < 5 Å).
#
include("render.jl")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const dir = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files", "lignin")

const lignin = (selection = "not resname TIP3 ACT", style = "Licorice 0.3 12 12", color = "Chain")
const acetone = (selection = "resname ACT", style = "VDW 1.0 12", color = "ColorID 7", material = "Transparent")
const water = (selection = "resname TIP3", style = "VDW 1.0 12", color = "ColorID 10", material = "Transparent")
cut(rep) = merge(rep, (selection = "($(rep.selection)) and x > 0",))
slab(rep) = merge(rep, (selection = "($(rep.selection)) and abs(z) < 5",))
# the lignin as spheres in the slabs
const lignin_vdw = merge(lignin, (style = "VDW 1.0 12",))

render(joinpath(dir, "particle_check.pdb"), joinpath(images, "lignin_check.webp");
    reps = [lignin], rotate = [(:z, 30), (:x, -70)])
render(joinpath(dir, "particle.pdb"), joinpath(images, "lignin_particle.webp");
    reps = [lignin], rotate = [(:z, 30), (:x, -70)])
render(joinpath(dir, "lignin.pdb"), joinpath(images, "lignin.webp");
    reps = [lignin, cut(acetone), cut(water)], rotate = [(:z, 30), (:x, -70)], cell = :origin)
render(joinpath(dir, "lignin.pdb"), joinpath(images, "lignin_slab.webp");
    reps = [slab(lignin_vdw), slab(acetone), slab(water)], cell = :origin)
render(joinpath(dir, "dispersed.pdb"), joinpath(images, "lignin_dispersed.webp");
    reps = [lignin, cut(acetone), cut(water)], rotate = [(:z, 30), (:x, -70)], cell = :origin)
render(joinpath(dir, "dispersed.pdb"), joinpath(images, "lignin_dispersed_slab.webp");
    reps = [slab(lignin_vdw), slab(acetone), slab(water)], cell = :origin)
