#
# Images of the minimal cell tutorial (docs/src/tutorials/cell.md), from the outputs of
# `cell.inp` and `solvation.inp`, in the directory where `prepare.jl` was run:
#
#     julia --project=dev/large_tests dev/tutorials/cell_images.jl path/to/dir
#
# The cell is cut as a cake (as cake.jl, with the view of `cake_rotate`): a quarter is
# removed, atom by atom. The outputs have up to 40 million beads, too many for VMD, and
# most of them cannot be seen: only the beads within `layer` Å of the two faces of the
# cut, and those within `surface_depth` Å of the outer surface of the membrane facing
# the viewer, are written
# to the PDB file that is rendered. The outputs are read line by line, not as a whole.
#
# The chain identifiers of that file mark the components: lipids (L) and proteins (P) of
# the membrane, the chromosome (D), the proteins (C), RNAs (R) and metabolites (S) of
# the cytosol, water (W) and glycerol (G).
#
include("cake.jl")      # cake_rotate, toward_viewer
include("render.jl")
using Printf: @printf
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const dir = ARGS[1]
const layer = 20.0
const surface_depth = 60.0   # the membrane is not a perfect sphere

# The center of the cell: that of the spheres of cell.inp
const center = let
    l = first(l for l in eachline(joinpath(dir, "cell.inp")) if occursin("inside sphere", l))
    parse.(Float64, split(l)[3:5])
end
count_atoms(file) = count(l -> startswith(l, "ATOM") || startswith(l, "HETATM"), eachline(file))
const n_membrane = count_atoms(joinpath(dir, "membrane.cif"))
const n_chromosome = count_atoms(joinpath(dir, "chromosome.cif"))
const lipids = Set(["CHOL", "CDL2", "DPSM", "POPC", "POPE", "POPG", "DPPC", "POPS", "DOPC", "DPPE", "PPCS"])
const rna = Set(["RA", "RU", "RG", "RC", "RA3", "RA5", "RU3", "RU5", "RG3", "RG5", "RC3", "RC5"])

# The component of an atom, from its index and names
function component(i, resname, name)
    i <= n_membrane && return resname in lipids ? 'L' : 'P'
    i <= n_membrane + n_chromosome && return 'D'
    resname == "W" && return 'W'
    resname == "GLY" && name in ("G1", "G2") && return 'G'
    resname in rna && return 'R'
    (name == "BB" || startswith(name, "SC")) && return 'C'
    return 'S'
end

# The fields of an atom line of the mmCIF outputs (atom names with spaces are quoted)
function fields(l)
    f = String[]
    i, n = 1, ncodeunits(l)
    while i <= n
        while i <= n && l[i] == ' '; i += 1; end
        i > n && break
        if l[i] == '\'' || l[i] == '"'
            q = l[i]; j = findnext(q, l, i + 1)
            push!(f, l[i+1:j-1]); i = j + 1
        else
            j = something(findnext(' ', l, i), n + 1)
            push!(f, l[i:j-1]); i = j
        end
    end
    return f
end

# The radius of the outer surface of the membrane: of its lipids (a few beads of the
# membrane proteins lie far out), a high percentile
const r_outer = let r = Float64[]
    for l in eachline(joinpath(dir, "membrane.cif"))
        (startswith(l, "ATOM") || startswith(l, "HETATM")) || continue
        f = fields(l)
        f[6] in lipids || continue
        push!(r, sqrt(sum(abs2, parse.(Float64, f[11:13]) .- center)))
    end
    sort!(r)[round(Int, 0.995 * length(r))]
end

function visible(x, y, z)
    in_cut(x, y) && return false
    x^2 + y^2 + z^2 > (r_outer + layer)^2 && return false   # stray beads
    (y < 0 && x > -layer) && return true          # the face x = 0 of the cut
    (x > 0 && y < layer) && return true           # the face y = 0 of the cut
    r = sqrt(x^2 + y^2 + z^2)
    return r > r_outer - surface_depth && x * toward_viewer[1] + y * toward_viewer[2] + z * toward_viewer[3] > 0
end

function visible_pdb(output)
    pdb = tempname() * ".pdb"
    n, kept = 0, 0
    open(pdb, "w") do io
        for l in eachline(output)
            (startswith(l, "ATOM") || startswith(l, "HETATM")) || continue
            n += 1
            f = fields(l)
            x, y, z = parse.(Float64, f[11:13]) .- center
            visible(x, y, z) || continue
            kept += 1
            c = component(n, f[6], f[4])
            @printf(io, "ATOM  %5d %-4s %-4s%1s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
                kept % 100000, first(f[4], 4), first(f[6], 4), c, 1, x, y, z)
        end
    end
    println("$output: $kept of $n atoms kept (outer radius $(round(r_outer; digits = 1)) A)")
    return pdb
end

reps = [
    (selection = "chain L", style = "VDW 2.4 8", color = "ColorID 3"),
    (selection = "chain P", style = "VDW 2.4 8", color = "ColorID 0"),
    (selection = "chain D", style = "VDW 2.4 8", color = "ColorID 4"),
    (selection = "chain C", style = "VDW 2.4 8", color = "ColorID 7"),
    (selection = "chain R", style = "VDW 2.4 8", color = "ColorID 1"),
    (selection = "chain S", style = "VDW 2.4 8", color = "ColorID 11"),
    (selection = "chain W", style = "VDW 2.4 8", color = "ColorID 15"),
    (selection = "chain G", style = "VDW 2.4 8", color = "ColorID 9"),
]
for (file, name) in (("cell.cif", "cell.webp"), ("cell_solvated.cif", "cell_solvated.webp"))
    pdb = visible_pdb(joinpath(dir, file))
    render(pdb, joinpath(images, name); reps, rotate = cake_rotate, size = (1600, 1600), width = 1000)
    rm(pdb)
end
