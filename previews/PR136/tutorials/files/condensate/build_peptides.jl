#
# Conformations of the elastin-like polypeptide (VPGVG)20 in the SIRAH coarse-grained
# representation (peptides.pdb, one model per conformation):
#
#   julia build_peptides.jl peptides.pdb
#
# Each chain grows residue by residue, as a self-avoiding random coil, with the dihedral
# angles drawn from the regions of the Ramachandran plot populated by disordered proteins
# (polyproline II, β and α). A residue is placed only if its atoms are at least 3.5 Å
# from those of the residues before its neighbor; at a dead end, the chain goes back by
# up to 5 residues, and after 200 dead ends it starts again. The atoms are then mapped to
# the beads of SIRAH, each placed on an atom: GN on N, GC on CA and GO on O (the
# backbone), BCB on CB (valine) and BCG on CG (proline). Glycine has no side-chain bead.
#
using LinearAlgebra
using Printf
using Random

const sequence = repeat("VPGVG", 20)
const nconformations = 10
const min_distance = 3.5

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

# Dihedral angles (φ, ψ) of a residue: the centers of the regions, with a spread of 15°.
# Proline has φ ≈ -65°; glycine can also take the left-handed helix; the residue
# before a proline is not in the α region.
function dihedrals(rng, residue, next)
    regions = [(-75.0, 145.0), (-120.0, 130.0), (-65.0, -40.0)]
    if residue == 'G'
        push!(regions, (75.0, 30.0))
    end
    if next == 'P'
        regions = regions[1:2]
    end
    φ, ψ = rand(rng, regions)
    if residue == 'P'
        φ = -65.0
    end
    return φ + 15 * randn(rng), ψ + 15 * randn(rng)
end

# The atoms of a residue: N, CA, C and, for valine and proline, CB or CG (the atom of
# the side-chain bead). O is computed when the next residue is placed (it depends on ψ).
struct Residue
    N::Vector{Float64}
    CA::Vector{Float64}
    C::Vector{Float64}
    side::Union{Nothing,Vector{Float64}}
end

function side_chain(rng, aa, N, CA, C)
    aa == 'G' && return nothing
    CB = place(C, N, CA, 1.53, 110.5, -120.7)
    aa == 'P' || return CB
    return place(N, CA, CB, 1.50, 104.0, rand(rng, (-30.0, 30.0)))
end

atoms(r::Residue) = isnothing(r.side) ? (r.N, r.CA, r.C) : (r.N, r.CA, r.C, r.side)

# The distance between two atoms (without the temporary vector of `norm(x - y)`: this is
# the inner loop of the growth of the chains)
distance(x, y) = sqrt((x[1] - y[1])^2 + (x[2] - y[2])^2 + (x[3] - y[3])^2)

# True if the atoms `xs` are at least `min_distance` from the atoms of the residues
# `1:last` and from the carbonyl oxygens `O[1:last]`
function fits(xs, chain, O, last)
    for x in xs
        for j in 1:last
            for y in atoms(chain[j])
                distance(x, y) < min_distance && return false
            end
            distance(x, O[j]) < min_distance && return false
        end
    end
    return true
end

# One conformation: the residues and the carbonyl oxygens, or `nothing` if the chain
# reached too many dead ends (it is then grown again from the start)
function try_grow_chain(rng; max_dead_ends = 200)
    n = length(sequence)
    N = [0.0, 0.0, 0.0]
    CA = [1.458, 0.0, 0.0]
    C = [2.009, 1.420, 0.0]
    chain = [Residue(N, CA, C, side_chain(rng, sequence[1], N, CA, C))]
    O = Vector{Float64}[]
    tries = 0
    dead_ends = 0
    while length(chain) < n
        i = length(chain)
        prev = chain[i]
        _, ψ = dihedrals(rng, sequence[i], sequence[i+1])
        φ, _ = dihedrals(rng, sequence[i+1], get(sequence, i + 2, ' '))
        Oi = place(prev.N, prev.CA, prev.C, 1.231, 120.5, ψ + 180.0)
        N = place(prev.N, prev.CA, prev.C, 1.329, 116.2, ψ)
        CA = place(prev.CA, prev.C, N, 1.458, 121.7, 180.0)
        C = place(prev.C, N, CA, 1.525, 111.2, φ)
        new = Residue(N, CA, C, side_chain(rng, sequence[i+1], N, CA, C))
        if fits(atoms(new), chain, O, i - 1) && fits((Oi,), chain, O, i - 2)
            push!(O, Oi)
            push!(chain, new)
            tries = 0
        elseif (tries += 1) > 100
            # a dead end: back by up to 5 residues
            (dead_ends += 1) > max_dead_ends && return nothing
            for _ in 1:min(5, length(chain) - 1)
                pop!(chain)
                pop!(O)
            end
            tries = 0
        end
    end
    # the oxygen of the last residue
    last = chain[end]
    for _ in 1:100
        _, ψ = dihedrals(rng, sequence[end], ' ')
        Oi = place(last.N, last.CA, last.C, 1.231, 120.5, ψ + 180.0)
        if fits((Oi,), chain, O, n - 2)
            push!(O, Oi)
            return chain, O
        end
    end
    return nothing
end

function grow_chain(rng)
    while true
        result = try_grow_chain(rng)
        isnothing(result) || return result
    end
end

const side_bead = Dict('V' => "BCB", 'P' => "BCG")

function write_peptides(file; seed = 1)
    rng = Xoshiro(seed)
    open(file, "w") do io
        println(io, "REMARK   (VPGVG)20, SIRAH coarse-grained, $nconformations conformations")
        for model in 1:nconformations
            chain, O = grow_chain(rng)
            println(io, @sprintf("MODEL     %4d", model))
            index = 0
            for (i, residue) in enumerate(chain)
                beads = [("GN", residue.N), ("GC", residue.CA)]
                if !isnothing(residue.side)
                    push!(beads, (side_bead[sequence[i]], residue.side))
                end
                push!(beads, ("GO", O[i]))
                for (name, x) in beads
                    index += 1
                    println(io, @sprintf("ATOM  %5d %-4s s%-2s A%4d    %8.3f%8.3f%8.3f  1.00  0.00",
                        index, name, sequence[i], i, x[1], x[2], x[3]))
                end
            end
            println(io, "ENDMDL")
            println("conformation $model")
        end
        println(io, "END")
    end
end

write_peptides(ARGS[1])
