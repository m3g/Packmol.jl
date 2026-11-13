#
# Images of the supercapacitor tutorial (docs/src/tutorials/supercapacitor.md): the
# initial placement (the `check_initial_placement` line of the input uncommented) and the
# packed system, seen along y, a slice of 12 Å in y, with whole ions. Run from the root of
# the repository, with the packmol executable built (juliac/build.jl):
#
#     julia --project=docs dev/tutorials/supercapacitor.jl
#
include("render.jl")
using PDBTools

const exe = joinpath(@__DIR__, "..", "..", "juliac", "build", "bin", "packmol")
const files = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files", "supercapacitor")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const work = mktempdir()
const ymax = 12.0

for f in readdir(files)
    cp(joinpath(files, f), joinpath(work, f); force = true)
end
input = read(joinpath(work, "supercapacitor.inp"), String)
write(joinpath(work, "check_supercapacitor.inp"), replace(input,
    "# check_initial_placement" => "check_initial_placement",
    "output supercapacitor.pdb" => "output check_supercapacitor.pdb",
))
for inp in ("check_supercapacitor.inp", "supercapacitor.inp")
    run(pipeline(`$exe -i $(joinpath(work, inp))`; stdout = devnull))
end

# The ions (whole) with the center in the slice 0 ≤ y < ymax, and the electrode atoms in
# it, written with the cell of the slice
function slice(file)
    atoms = read_pdb(file)
    keep = Atom[]
    for mol in eachresidue(atoms)
        y = sum(at.y for at in mol) / length(mol)
        if resname(mol) in ("CAT", "ANO")
            append!(keep, [at for at in mol if 0 <= at.y < ymax])
        elseif 0 <= y < ymax
            append!(keep, collect(mol))
        end
    end
    s = Structure(keep)
    s.unitcell = UnitCell(39.36, ymax, 100.0, 90.0, 90.0, 90.0)
    out = tempname() * ".pdb"
    write_pdb(out, s)
    return out
end

reps = [
    (selection = "resname CAT ANO", style = "VDW 0.7 12", color = "ColorID 2"),
    (selection = "resname BMIM", style = "VDW 1.0 12", color = "ColorID 0"),
    (selection = "resname PF6", style = "VDW 1.0 12", color = "ColorID 1"),
]
for (file, image) in (("check_supercapacitor.pdb", "supercapacitor_check.webp"), ("supercapacitor.pdb", "supercapacitor.webp"))
    render(slice(joinpath(work, file)), joinpath(images, image);
        reps, rotate = [(:x, -90), (:z, -90)], cell = :unitcell, size = (1800, 1000), width = 900)
end
