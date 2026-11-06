#
# Images of the exosome tutorial (docs/src/tutorials/exosome.md), from the outputs of
# `exosome.inp` of the packmol_exosome set of MolSimTestFiles (with and without
# `check_initial_placement`), in the directory of the set, after `proteins.inp`:
#
#     julia --project=dev/large_tests dev/tutorials/exosome_images.jl exosome.inp exosome_check.cif exosome.cif
#
# The systems have 13 million atoms: the molecules with an atom in a slab through the
# center of the vesicle (|y| < 10 Å) are written first to a smaller PDB file, which is
# rendered. The chain identifiers of that file mark the membrane proteins (M) and the
# contents of the lumen (N) of the fixed `proteins.cif` (the output of `proteins.inp`).
#
using PDBToolsBase: read_mmcif, write_pdb, Atom
include("render.jl")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const input, outputs = ARGS[1], ARGS[2:end]
const slab = 10.0

# The structure files and numbers of molecules of an input
function structures(input)
    types = Tuple{String,Int}[]
    file = ""
    for line in eachline(input)
        w = split(strip(first(split(line, '#'))))
        isempty(w) && continue
        w[1] == "structure" && (file = String(w[2]))
        w[1] == "number" && push!(types, (file, parse(Int, w[2])))
    end
    return types
end
types = structures(input)
natoms = Dict(f => length(read_mmcif(joinpath(dirname(input), f))) for f in unique(first.(types)))
# The atoms of the membrane proteins, the first ones of proteins.cif
n_membrane = sum(n * length(read_mmcif(joinpath(dirname(input), f)))
    for (f, n) in structures(joinpath(dirname(input), "proteins.inp")) if f in ("cd81.cif", "cd9.cif"))

function slab_file(output)
    atoms = read_mmcif(output)
    kept = Atom[]
    i = 0
    for (file, n) in types
        if file == "proteins.cif"
            # the fixed proteins: their atoms in the slab
            for (k, a) in enumerate(@view(atoms[i+1:i+natoms[file]]))
                abs(a.y) < slab || continue
                b = copy(a); b.chain = k <= n_membrane ? "M" : "N"
                push!(kept, b)
            end
            i += natoms[file]
            continue
        end
        for _ in 1:n
            mol = @view(atoms[i+1:i+natoms[file]])
            any(a -> abs(a.y) < slab, mol) && append!(kept, mol)
            i += natoms[file]
        end
    end
    i == length(atoms) || error("$output has $(length(atoms)) atoms, the input $i")
    pdb = tempname() * ".pdb"
    write_pdb(pdb, kept)
    return pdb
end

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
    pdb = slab_file(file)
    render(pdb, joinpath(images, name); reps, rotate = [(:x, -90)])
    rm(pdb)
end
