#
# Images of the lipid nanoparticle tutorial (docs/src/tutorials/lnp.md), from the outputs
# of the input files of the packmol_lnp set of MolSimTestFiles, in the directory of the
# set:
#
#     julia --project=dev/large_tests dev/tutorials/lnp_images.jl monolayers.inp monolayers.cif lnp_monolayers.webp
#     julia --project=dev/large_tests dev/tutorials/lnp_images.jl lnp.inp lnp_check.cif lnp_check.webp
#     julia --project=dev/large_tests dev/tutorials/lnp_images.jl lnp.inp lnp.cif lnp.webp
#
# The nanoparticle is cut as a cake (cake.jl): the molecules with the center in the
# quarter x > 0, y < 0 are removed (the RNA, a single molecule each, bead by bead), and
# the water outside the particle is kept behind it. The molecules kept are written first
# to a smaller PDB file, which is rendered, with the chain identifiers marking the
# components: M the lipids of the monolayers, F the POPC of the lipid phase, R the RNA,
# C the water of the compartments, O the water outside.
#
include("cake.jl")
include("render.jl")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const input, output, image = ARGS

# The component of each structure block, and how it is cut
function tag(block)
    block.atoms && return ("M", :whole)
    startswith(block.file, "rna") && return ("R", :split)
    block.file == "wt4.pdb" && return outside_radius(block) > 0 ? ("O", :outside) : ("C", :whole)
    return ("F", :whole)
end

heads = "name BCO BPE BFO1"
reps = [
    (selection = "chain M and $heads", style = "VDW 2.4 8", color = "ColorID 0"),
    (selection = "chain M and not $heads", style = "VDW 1.8 6", color = "ColorID 3"),
    (selection = "chain F", style = "VDW 1.8 6", color = "ColorID 5"),
    (selection = "chain R", style = "VDW 2.4 8", color = "ColorID 1"),
    (selection = "chain C", style = "VDW 1.8 6", color = "ColorID 15"),
    (selection = "chain O", style = "VDW 1.8 6", color = "ColorID 15"),
]
pdb = cake_file(input, output; tag, ahead = 70.0)
render(pdb, joinpath(images, image); reps, rotate = cake_rotate, size = (1600, 1600), width = 1000)
rm(pdb)
