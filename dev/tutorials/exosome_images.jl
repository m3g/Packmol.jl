#
# Images of the exosome tutorial (docs/src/tutorials/exosome.md), from the outputs of
# `exosome.inp` of the packmol_exosome set of MolSimTestFiles (with and without
# `check_initial_placement`), in the directory of the set, after `proteins.inp`:
#
#     julia --project=dev/large_tests dev/tutorials/exosome_images.jl exosome.inp exosome_check.cif exosome.cif
#
# The vesicle is cut as a cake (cake.jl): a quarter of it is removed (the fixed proteins
# atom by atom), and the solvent outside the vesicle is kept behind it. The heavy atoms
# kept (of the 13 million atoms) are written first to a smaller PDB file, which is
# rendered. The chain identifiers of that file mark the membrane proteins (M) and the
# contents of the lumen (N) of the fixed `proteins.cif` (the output of `proteins.inp`),
# and the solvent outside the vesicle (O).
#
include("cake.jl")
include("render.jl")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const input, outputs = ARGS[1], ARGS[2:end]

# The atoms of the membrane proteins, the first ones of proteins.cif
const n_membrane = sum(b.number * length(read_mmcif(joinpath(dirname(input), b.file)))
    for b in structure_blocks(joinpath(dirname(input), "proteins.inp")) if b.file in ("cd81.cif", "cd9.cif"))

# The component of each structure block, and how it is cut: the solvent outside the
# vesicle is the one outside the sphere of 258.1 Å
function tag(block)
    block.file == "proteins.cif" && return (k -> k <= n_membrane ? "M" : "N", :split)
    return outside_radius(block) > 250 ? ("O", :outside) : ("X", :whole)
end
heavy(a) = !startswith(a.name, "H")

nucleic = "resname A U G C DA DT DG DC"
reps = [
    (selection = "chain M and noh", style = "VDW 1.0 12", color = "ColorID 0"),
    (selection = "chain N and not ($nucleic) and noh", style = "VDW 1.0 12", color = "ColorID 7"),
    (selection = "chain N and $nucleic and noh", style = "VDW 1.0 12", color = "ColorID 4"),
    (selection = "resname POPC POPE POPS DPPC and noh", style = "VDW 1.0 8", color = "ColorID 3"),
    (selection = "resname TIP3 and noh", style = "VDW 1.0 8", color = "ColorID 15"),
    (selection = "resname SOD", style = "VDW 1.6 12", color = "ColorID 1"),
    (selection = "resname CLA", style = "VDW 1.6 12", color = "ColorID 9"),
]
for (file, name) in zip(outputs, ("exosome_check.webp", "exosome.webp"))
    pdb = cake_file(input, file; tag, ahead = 30.0, keep_atom = heavy)
    render(pdb, joinpath(images, name); reps, rotate = cake_rotate, size = (1600, 1600), width = 1000)
    rm(pdb)
end
