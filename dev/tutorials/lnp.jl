#
# The files of the lipid nanoparticle tutorial (docs/src/tutorials/lnp.md), the
# `packmol_lnp` set of MolSimTestFiles:
#
#     julia --project=dev/large_tests dev/tutorials/lnp.jl <directory>
#
# writes to <directory> the input files (docs/src/tutorials/files/lnp), the SIRAH lipids
# and water of the `packmol_zika` set (popc.pdb, pope.pdb, pops.pdb, wt4.pdb), and the
# pieces of RNA, rna1.cif, rna2.cif and rna3.cif: compact self-avoiding chains of
# nucleotides of 3 beads (BB1, BB2, SC1) in spheres, centered at the origin. They are
# placeholders, with the size of the RNA a compartment holds, not structures of a
# coarse-grained model of RNA.
#
using PDBToolsBase: Atom, write_mmcif
using LazyArtifacts: @artifact_str
using StaticArrays: SVector
using LinearAlgebra: normalize, cross, norm
using Random: Xoshiro

const dir = mkpath(ARGS[1])
const rng = Xoshiro(3)

# A chain of `nnt` nucleotides within a sphere of radius `radius`: the backbone beads
# `step` apart, all beads at least `rmin` apart; at a dead end, the chain goes back by
# up to 5 nucleotides
function rna_globule(nnt, radius; step = 6.0, rmin = 4.0)
    beads = SVector{3,Float64}[]
    path = SVector{3,Float64}[]
    tries = 0
    while length(beads) < 3nnt
        u = normalize(SVector(randn(rng), randn(rng), randn(rng)))
        q = isempty(path) ? zero(SVector{3,Float64}) : path[end] + step * u
        v = normalize(SVector(randn(rng), randn(rng), randn(rng)))
        new = [q, q + 3.5 * normalize(cross(u, v)), q + 4.5 * normalize(cross(u, v)) + 2.5 * u]
        if all(norm(x) <= radius for x in new) && all(norm(x - y) >= rmin for x in new for y in beads)
            append!(beads, new)
            push!(path, q)
            tries = 0
        elseif (tries += 1) > 200
            for _ in 1:min(5, length(path) - 1)
                pop!(path)
                resize!(beads, length(beads) - 3)
            end
            tries = 0
        end
    end
    c = sum(beads) / length(beads)
    names = repeat(["BB1", "BB2", "SC1"], nnt)
    bases = rand(rng, ["A", "U", "G", "C"], nnt)
    return [Atom(name = names[i], resname = bases[cld(i, 3)], chain = "R", resnum = cld(i, 3),
        x = beads[i][1] - c[1], y = beads[i][2] - c[2], z = beads[i][3] - c[3]) for i in eachindex(beads)]
end

# Nucleotides and radius of the sphere (Å) of the RNA of each compartment
for (i, (nnt, r)) in enumerate(((1500, 62.0), (1000, 54.0), (700, 47.0)))
    write_mmcif(joinpath(dir, "rna$i.cif"), rna_globule(nnt, r))
    println("rna$i.cif: $nnt nucleotides, radius $r Å")
end

# The packmol_zika set, of the Artifacts.toml of dev/large_tests (the macro looks for the
# file next to the source it is evaluated in)
const zika = include_string(@__MODULE__, "artifact\"packmol_zika\"",
    joinpath(@__DIR__, "..", "large_tests", "run.jl"))
for file in ("popc.pdb", "pope.pdb", "pops.pdb", "wt4.pdb")
    cp(joinpath(zika, file), joinpath(dir, file); force = true)
    chmod(joinpath(dir, file), 0o644)
end
const inputs = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files", "lnp")
for file in readdir(inputs)
    cp(joinpath(inputs, file), joinpath(dir, file); force = true)
end
