#
# Images of the gold nanoparticle tutorial (docs/src/tutorials/gold_nanoparticle.md): the
# initial placement (the `check_initial_placement` line of the input uncommented) and the
# packed system, a slab of 16 Å through the center (|z| < 8, whole molecules, and the gold
# atoms in it), seen along z. Run from the root of the repository, with the packmol
# executable built (juliac/build.jl):
#
#     julia --project=docs dev/tutorials/gold_nanoparticle.jl
#
include("render.jl")
using PDBTools

const exe = joinpath(@__DIR__, "..", "..", "juliac", "build", "bin", "packmol")
const files = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files", "gold_nanoparticle")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const work = mktempdir()
const half = 8.0

for f in readdir(files)
    cp(joinpath(files, f), joinpath(work, f); force = true)
end
input = read(joinpath(work, "gold_nanoparticle.inp"), String)
write(joinpath(work, "check_gold_nanoparticle.inp"), replace(input,
    "# check_initial_placement" => "check_initial_placement",
    "output gold_nanoparticle.pdb" => "output check_gold_nanoparticle.pdb",
))
for inp in ("check_gold_nanoparticle.inp", "gold_nanoparticle.inp")
    run(pipeline(`$exe -i $(joinpath(work, inp))`; stdout = devnull))
end

# The molecules (whole) with the center in the slab |z| < half, and the gold atoms in it
function slab(file)
    atoms = read_pdb(file)
    keep = Atom[]
    for mol in eachresidue(atoms)
        if resname(mol) == "AU"
            append!(keep, [at for at in mol if abs(at.z) < half])
        elseif abs(sum(at.z for at in mol) / length(mol)) < half
            append!(keep, collect(mol))
        end
    end
    out = tempname() * ".pdb"
    write_pdb(out, keep)
    return out
end

reps = [
    (selection = "resname AU", style = "VDW 1.0 20", color = "ColorID 4"),
    (selection = "resname DDT and element S", style = "VDW 1.0 16", color = "ColorID 3"),
    (selection = "resname DDT and not element S", style = "Licorice 0.35 12 12", color = "ColorID 2"),
    (selection = "resname TOL", style = "Licorice 0.25 12 12", color = "ColorID 10", material = "Transparent"),
]
for (file, image) in (("check_gold_nanoparticle.pdb", "gold_nanoparticle_check.webp"), ("gold_nanoparticle.pdb", "gold_nanoparticle.webp"))
    render(slab(joinpath(work, file)), joinpath(images, image); reps, size = (1600, 1600), width = 800)
end
