#
# The "cake" views of the spherical particles of the tutorials (the lipid nanoparticle,
# the Zika virion, the exosome): a quarter of the particle is cut out (the molecules with
# the center in x > 0, y < 0; fixed molecules and other large ones atom by atom), and the
# solvent outside the particle is kept only behind it, cut by a plane turned by `turn`
# degrees from the direction of the viewer, `ahead` Å in front of the center. Thus the
# inside of the particle is seen in the cut, its spherical shape in front, and the
# periodic cell of solvent behind it.
#
# `cake_file(input, output; tag, ahead)` writes the atoms kept to a temporary PDB file,
# rendered with `render(pdb, image; reps, rotate = cake_rotate, ...)` (render.jl). `tag`
# is a function of a structure block of the input (see `structure_blocks`), returning
# the chain identifier of its atoms in that file (a string, a function of the index of
# the atom in the molecule, or `nothing` to keep the chains of the output) and how the
# block is cut: `:split` (atom by atom), `:whole` (by the center of the molecules) or
# `:outside` (by the center, and only behind the plane: the solvent outside the
# particle; by the plane only, not by the quarter, beyond the radius `uncut_beyond`).
#
using PDBToolsBase: read_mmcif, read_pdb, write_pdb, Atom

read_structure(file) = endswith(file, ".cif") ? read_mmcif(file) : read_pdb(file)

# The structure blocks of an input file: the structure file, the number of molecules,
# the words of the lines of the block that apply to whole molecules, and whether the
# block has `atoms` sections (constraints on some atoms)
function structure_blocks(input)
    blocks = @NamedTuple{file::String, number::Int, lines::Vector{Vector{String}}, atoms::Bool}[]
    file, number, lines, atoms, in_atoms = "", 0, Vector{String}[], false, false
    for line in eachline(input)
        w = String.(split(strip(first(split(line, '#')))))
        isempty(w) && continue
        if w[1] == "structure"
            file, number, lines, atoms = w[2], 0, Vector{String}[], false
        elseif w[1] == "number"
            number = parse(Int, w[2])
        elseif w[1] == "atoms"
            atoms = in_atoms = true
        elseif w[1] == "end" && length(w) > 1 && w[2] == "atoms"
            in_atoms = false
        elseif w[1] == "end" && length(w) > 1 && w[2] == "structure"
            push!(blocks, (; file, number, lines, atoms))
        elseif !in_atoms && !isempty(file)
            push!(lines, w)
        end
    end
    return blocks
end

# The radius of an `outside sphere` constraint of a block (0 if none)
function outside_radius(block)
    i = findfirst(w -> length(w) == 6 && w[1] == "outside" && w[2] == "sphere", block.lines)
    return isnothing(i) ? 0.0 : parse(Float64, block.lines[i][6])
end

# The view of the images and the plane cutting the solvent outside
const cake_rotate = [(:x, -90), (:y, -25), (:x, 35)]
# The direction towards the viewer and the right of the image, in the coordinates of
# the system, for the view of `cake_rotate`
const toward_viewer = (0.3462, -0.7424, 0.5736)
const image_right = (0.9063, 0.4226, 0.0)
in_cut(x, y) = x > 0 && y < 0
function behind_plane(x, y, z; turn, ahead)
    n = cosd(turn) .* toward_viewer .+ sind(turn) .* image_right
    return n[1] * x + n[2] * y + n[3] * z < ahead
end

function cake_file(input, output; tag, ahead, turn = 40.0, uncut_beyond = Inf, keep_atom = a -> true)
    blocks = structure_blocks(input)
    natoms = Dict(b.file => length(read_structure(joinpath(dirname(input), b.file))) for b in blocks)
    atoms = read_structure(output)
    kept = Atom[]
    i = 0
    for block in blocks
        chain, kind = tag(block)
        n = natoms[block.file]
        for _ in 1:block.number
            mol = @view(atoms[i+1:i+n])
            i += n
            if kind == :split
                for (k, a) in enumerate(mol)
                    (in_cut(a.x, a.y) || !keep_atom(a)) && continue
                    b = copy(a)
                    isnothing(chain) || (b.chain = chain isa String ? chain : chain(k))
                    push!(kept, b)
                end
            else
                xc = sum(a.x for a in mol) / n
                yc = sum(a.y for a in mol) / n
                zc = sum(a.z for a in mol) / n
                in_cut(xc, yc) && (kind != :outside || xc^2 + yc^2 + zc^2 < uncut_beyond^2) && continue
                kind == :outside && !behind_plane(xc, yc, zc; turn, ahead) && continue
                for a in mol
                    keep_atom(a) || continue
                    b = copy(a)
                    isnothing(chain) || (b.chain = chain)
                    push!(kept, b)
                end
            end
        end
    end
    i == length(atoms) || error("$output has $(length(atoms)) atoms, the input $i")
    pdb = tempname() * ".pdb"
    write_pdb(pdb, kept)
    return pdb
end
