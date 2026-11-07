#
# Images of the lipid nanoparticle tutorial (docs/src/tutorials/lnp.md), from the outputs
# of the input files of the packmol_lnp set of MolSimTestFiles, in the directory of the
# set:
#
#     julia --project=dev/large_tests dev/tutorials/lnp_images.jl monolayers.inp monolayers.cif lnp_monolayers.webp
#     julia --project=dev/large_tests dev/tutorials/lnp_images.jl lnp.inp lnp_check.cif lnp_check.webp
#     julia --project=dev/large_tests dev/tutorials/lnp_images.jl lnp.inp lnp.cif lnp.webp
#
# The nanoparticle is cut as a cake: the molecules with the center in the quarter x > 0,
# y < 0 are removed (the RNA, a single molecule each, bead by bead). Of the water outside
# the particle, only two walls of `wall` Å along the faces of the cut are kept, so that
# the surface of the particle is seen. The molecules kept are written first to a smaller
# PDB file, which is rendered, with the chain identifiers marking the components: M the
# lipids of the monolayers, F the POPC of the lipid phase, R the RNA, C the water of the
# compartments, O the water outside.
#
using PDBToolsBase: read_mmcif, read_pdb, write_pdb, Atom
include("render.jl")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const input, output, image = ARGS

# The structure file, number of molecules and component of each structure block
function structures(input)
    types = Tuple{String,Int,String}[]
    file, number, oriented, outside = "", 0, false, false
    for line in eachline(input)
        w = split(strip(first(split(line, '#'))))
        isempty(w) && continue
        if w[1] == "structure"
            file, oriented, outside = String(w[2]), false, false
        elseif w[1] == "number"
            number = parse(Int, w[2])
        elseif w[1] == "atoms"
            oriented = true
        elseif w[1] == "outside" && file == "wt4.pdb"
            outside = true
        elseif w[1] == "end" && length(w) > 1 && w[2] == "structure"
            tag = oriented ? "M" : startswith(file, "rna") ? "R" :
                  file == "wt4.pdb" ? (outside ? "O" : "C") : "F"
            push!(types, (file, number, tag))
        end
    end
    return types
end

read_structure(file) = endswith(file, ".cif") ? read_mmcif(file) : read_pdb(file)
cut(x, y) = x > 0 && y < 0 # the quarter removed
const wall = 25.0
wall_of_water(x, y) = (x > 0 && 0 <= y < wall) || (y < 0 && -wall < x <= 0)

function cake_file(input, output)
    types = structures(input)
    natoms = Dict(f => length(read_structure(joinpath(dirname(input), f))) for f in unique(first.(types)))
    atoms = read_mmcif(output)
    kept = Atom[]
    i = 0
    for (file, n, tag) in types
        for _ in 1:n
            mol = @view(atoms[i+1:i+natoms[file]])
            i += natoms[file]
            if tag == "R"
                for a in mol
                    cut(a.x, a.y) && continue
                    b = copy(a); b.chain = tag; push!(kept, b)
                end
            else
                xc = sum(a.x for a in mol) / length(mol)
                yc = sum(a.y for a in mol) / length(mol)
                cut(xc, yc) && continue
                tag == "O" && !wall_of_water(xc, yc) && continue
                for a in mol
                    b = copy(a); b.chain = tag; push!(kept, b)
                end
            end
        end
    end
    i == length(atoms) || error("$output has $(length(atoms)) atoms, the input $i")
    pdb = tempname() * ".pdb"
    write_pdb(pdb, kept)
    return pdb
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
pdb = cake_file(input, output)
render(pdb, joinpath(images, image); reps,
    rotate = [(:x, -90), (:y, -25), (:x, 35)], size = (1600, 1600), width = 1000)
rm(pdb)
