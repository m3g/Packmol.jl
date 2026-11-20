#
# Conformations of a thiol-terminated polystyrene, sec-butyl-[CH2-CH(C6H5)]20-CH2CH2-S⁻
# (anionic polystyrene end-capped with ethylene sulfide), all atoms, written to a PDB
# file with one model per conformation:
#
#   julia build_polystyrene.jl polystyrene.pdb
#
# A proof of concept only, not a force field: the chain grows from the sulfur atom along
# its backbone, with the torsions drawn from the trans and gauche states (180°, 60°,
# -60°, with a spread of 10°), atactic (each phenyl group on a random side of the
# backbone), as a self-avoiding walk: an atom is placed only if it is at least 2.2 Å from
# the atoms of the chain more than three bonds away from it. The geometry: C-C 1.53 Å, C-S
# 1.82 Å, C-H 1.09 Å, aromatic C-C 1.39 Å, C-H 1.08 Å, tetrahedral angles of 109.5°
# (112° along the backbone), the phenyl rings perpendicular to the plane of their
# backbone carbon atoms.
#
# The atoms are in the order of the chain from the sulfur atom, in the groups that are
# the monomers of the input file: S-CH2 (4 atoms), CH2 (3), then 20 times CH(C6H5) (13)
# and CH2 (3), and the sec-butyl group (13). Each group starts with its backbone carbon
# atom (S in the first one), followed by its substituents.
#
using LinearAlgebra
using Printf
using Random

const nunits = 20
const nconformations = 10
const min_distance = 2.2

# The position of an atom bonded to `c`, with the bond length, the angle b-c-new and
# the dihedral a-b-c-new (angles in degrees)
function place(a, b, c, bond, angle, dihedral)
    bc = normalize(c - b)
    n = normalize(cross(b - a, bc))
    m = cross(n, bc)
    θ = deg2rad(angle)
    φ = deg2rad(dihedral)
    return c + bond * (-cos(θ) * bc + sin(θ) * cos(φ) * m + sin(θ) * sin(φ) * n)
end

# The two tetrahedral directions of the substituents of the atom at `c`, bonded to `p`
# and `n` along the backbone: in the plane bisecting p-c-n, on either side of it
function substituent_directions(p, c, n)
    u = normalize(p - c)
    v = normalize(n - c)
    bisector = -normalize(u + v)
    normal = normalize(cross(u, v))
    θ = deg2rad(109.5 / 2)
    return bisector * cos(θ) + normal * sin(θ), bisector * cos(θ) - normal * sin(θ)
end

# A torsion of the backbone: trans or gauche, with a spread of 10°
torsion(rng) = rand(rng, (180.0, 180.0, 60.0, -60.0)) + 10 * randn(rng)

# The phenyl group bonded to the atom at `c` along `direction`: the ring carbon atoms
# (ipso, ortho, meta, para, meta, ortho) and the hydrogen atoms (ortho to ortho), the
# ring perpendicular to the plane of `p`, `c` and `n`, and to the backbone
function phenyl(p, c, n, direction)
    ipso = c + 1.51 * direction
    center = ipso + 1.39 * direction
    # the ring in the plane of `direction` and of the normal to the backbone plane
    normal = normalize(cross(p - c, n - c))
    w = normalize(normal - dot(normal, direction) * direction)
    ring(k) = -cos(deg2rad(60k)) * direction + sin(deg2rad(60k)) * w
    carbons = [center + 1.39 * ring(k) for k in 0:5]
    hydrogens = [center + (1.39 + 1.08) * ring(k) for k in 1:5]
    return carbons, hydrogens
end

# The backbone: S, the two CH2 of the end group, the 2 × 20 carbon atoms of the units,
# and the three carbon atoms of the sec-butyl group along its ethyl branch (CH, CH2, CH3)
const nbackbone = 3 + 2 * nunits + 3

# An atom: its name, residue name, residue number, element, position, and the backbone
# atom it hangs from with the number of bonds to it (to count the bonds between atoms)
struct ChainAtom
    name::String
    resname::String
    resnum::Int
    element::String
    x::Vector{Float64}
    anchor::Int
    depth::Int
end

# The residue of the backbone atom `c`: the end group (1), the styrene units (2 to 21),
# the sec-butyl group (22)
residue(c) = c <= 3 ? ("ETS", 1) : c <= nbackbone - 3 ? ("STY", (c - 2) ÷ 2 + 1) : ("SBU", nunits + 2)

# The substituents of the backbone atom `c`, bonded to `p` and `n` along the backbone (`n`
# is `nothing` for the last one)
function substituents(rng, c, x, p, n, pp)
    resname, resnum = residue(c)
    atoms = ChainAtom[]
    add!(name, element, r, depth) = push!(atoms, ChainAtom(name, resname, resnum, element, r, c, depth))
    if c == 1
        return atoms  # the sulfur atom: bound to the gold
    elseif c in (2, 3)
        d1, d2 = substituent_directions(p, x, n)
        add!("H$(c-1)A", "H", x + 1.09 * d1, 1)
        add!("H$(c-1)B", "H", x + 1.09 * d2, 1)
    elseif c <= nbackbone - 3 && iseven(c)
        # CH(C6H5): atactic, the phenyl group on either side of the backbone
        d1, d2 = substituent_directions(p, x, n)
        dphenyl, dh = rand(rng, Bool) ? (d1, d2) : (d2, d1)
        add!("HA", "H", x + 1.09 * dh, 1)
        carbons, hydrogens = phenyl(p, x, n, dphenyl)
        for (name, r, depth) in zip(("CG", "CD1", "CE1", "CZ", "CE2", "CD2"), carbons, (1, 2, 3, 4, 3, 2))
            add!(name, "C", r, depth)
        end
        for (name, r, depth) in zip(("HD1", "HE1", "HZ", "HE2", "HD2"), hydrogens, (3, 4, 5, 4, 3))
            add!(name, "H", r, depth)
        end
    elseif c <= nbackbone - 3
        # CH2
        d1, d2 = substituent_directions(p, x, n)
        add!("HB1", "H", x + 1.09 * d1, 1)
        add!("HB2", "H", x + 1.09 * d2, 1)
    elseif c == nbackbone - 2
        # the CH of the sec-butyl group, with its methyl group
        d1, d2 = substituent_directions(p, x, n)
        dmethyl, dh = rand(rng, Bool) ? (d1, d2) : (d2, d1)
        add!("H1", "H", x + 1.09 * dh, 1)
        methyl = x + 1.53 * dmethyl
        add!("C2", "C", methyl, 1)
        for (k, φ) in enumerate((60.0, 180.0, 300.0))
            add!("H2$k", "H", place(p, x, methyl, 1.09, 109.5, φ), 2)
        end
    elseif c == nbackbone - 1
        d1, d2 = substituent_directions(p, x, n)
        add!("H3A", "H", x + 1.09 * d1, 1)
        add!("H3B", "H", x + 1.09 * d2, 1)
    else
        # the terminal CH3
        for (k, φ) in enumerate((60.0, 180.0, 300.0))
            add!("H4$k", "H", place(pp, p, x, 1.09, 109.5, φ), 1)
        end
    end
    return atoms
end

# The name of the backbone atom `c`
function backbone_atom(c, x)
    resname, resnum = residue(c)
    name, element = if c == 1
        "S1", "S"
    elseif c <= 3
        "C$(c-1)", "C"
    elseif c <= nbackbone - 3
        iseven(c) ? "CA" : "CB", "C"
    elseif c == nbackbone - 2
        "C1", "C"  # the CH of the sec-butyl group (its methyl group is C2)
    elseif c == nbackbone - 1
        "C3", "C"
    else
        "C4", "C"
    end
    return ChainAtom(name, resname, resnum, element, x, c, 0)
end

# The number of bonds between two atoms (along the backbone and the substituents; the
# atoms hanging from the same backbone atom are in the same rigid group)
bonds_between(a, b) = abs(a.anchor - b.anchor) + a.depth + b.depth

# Whether the atoms `new` are at least `min_distance` from the atoms of `chain` more than
# three bonds away
function fits(new, chain)
    for a in new, b in chain
        a.anchor == b.anchor && continue
        bonds_between(a, b) > 3 || continue
        norm(a.x - b.x) < min_distance && return false
    end
    return true
end

# One conformation: the atoms, grouped by backbone atom (the backbone atom first, then its
# substituents), or `nothing` if the walk reached too many dead ends
function try_conformation(rng; max_dead_ends = 200)
    x = [[0.0, 0.0, 0.0], [1.82, 0.0, 0.0], [2.43, 1.72, 0.0]]
    groups = [[backbone_atom(c, x[c])] for c in 1:3]
    append!(groups[2], substituents(rng, 2, x[2], x[1], x[3], nothing))
    dead_ends = 0
    tries = 0
    while length(x) < nbackbone
        c = length(x) + 1
        new_x = place(x[c-3], x[c-2], x[c-1], 1.53, 112.0, torsion(rng))
        new_group = [backbone_atom(c, new_x)]
        # the substituents of the previous backbone atom, now that its neighbors are known
        previous = substituents(rng, c - 1, x[c-1], x[c-2], new_x, x[c-3])
        placed = reduce(vcat, groups)
        if fits(vcat(new_group, previous), placed) && fits(new_group, previous)
            append!(groups[c-1], previous)
            push!(x, new_x)
            push!(groups, new_group)
            tries = 0
        elseif (tries += 1) > 30
            # a dead end: back by up to 4 backbone atoms
            (dead_ends += 1) > max_dead_ends && return nothing
            for _ in 1:min(4, length(x) - 3)
                pop!(x)
                pop!(groups)
            end
            # the substituents of the new last atom are placed again with the next one
            groups[end] = groups[end][1:1]
            tries = 0
        end
    end
    # the hydrogen atoms of the terminal methyl group
    last = substituents(rng, nbackbone, x[end], x[end-1], nothing, x[end-2])
    fits(last, reduce(vcat, groups)) || return nothing
    append!(groups[end], last)
    return reduce(vcat, groups)
end

function conformation(rng)
    while true
        atoms = try_conformation(rng)
        isnothing(atoms) || return atoms
    end
end

function write_polystyrene(file; seed = 1)
    rng = Xoshiro(seed)
    open(file, "w") do io
        println(io, "REMARK   sec-butyl-(styrene)$nunits-CH2CH2-S-, $nconformations conformations")
        for model in 1:nconformations
            atoms = conformation(rng)
            println(io, @sprintf("MODEL     %4d", model))
            for (index, a) in enumerate(atoms)
                println(io, @sprintf("ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s",
                    index, a.name, a.resname, a.resnum, a.x[1], a.x[2], a.x[3], a.element))
            end
            println(io, "ENDMDL")
            println("conformation $model: $(length(atoms)) atoms")
        end
        println(io, "END")
    end
end

write_polystyrene(ARGS[1])
