#
# Images of the Zika virion tutorial (docs/src/tutorials/zika.md), from the outputs of
# the input of the packmol_zika set of MolSimTestFiles (with and without
# `check_initial_placement`), in the directory of the set:
#
#     julia --project=dev/large_tests dev/tutorials/zika.jl zika.inp zika_check.pdb zika.pdb
#
# The virion is cut as a cake (cake.jl): a quarter of it is removed (the capsid, fixed,
# atom by atom), and the solvent outside the membrane is kept behind it. The molecules
# kept are written first to a smaller PDB file, which is rendered, with the chain
# identifiers marking the components: those of the capsid (proteins E and M, colored by
# chain), B the lipids, W the solvent inside, Z the solvent outside.
#
include("cake.jl")
include("render.jl")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const input, outputs = ARGS[1], ARGS[2:end]

# The component of each structure block, and how it is cut: the solvent outside the
# membrane is the one outside spheres of 170 Å or more
function tag(block)
    block.file == "capsid.pdb" && return (nothing, :split)
    startswith(block.file, "pop") && return ("B", :whole)
    return outside_radius(block) >= 170 ? ("Z", :outside) : ("W", :whole)
end

# The chains of the capsid: K, M, O the protein E, L, N, P the protein M
chain_colors = ("K" => 0, "M" => 7, "O" => 4, "L" => 11, "N" => 9, "P" => 10)
reps = [
    [(selection = "chain $c", style = "VDW 1.6 12", color = "ColorID $id") for (c, id) in chain_colors]...,
    (selection = "chain B", style = "VDW 1.4 12", color = "ColorID 3"),
    (selection = "chain W Z and resname WT4 WLS", style = "VDW 1.0 8", color = "ColorID 15"),
    (selection = "resname NaW ClW", style = "VDW 1.6 12", color = "ResName"),
]
for (file, name) in zip(outputs, ("zika_check.webp", "zika.webp"))
    pdb = cake_file(input, file; tag, ahead = 0.0)
    render(pdb, joinpath(images, name); reps, rotate = cake_rotate, size = (1600, 1600), width = 1000)
    rm(pdb)
end
