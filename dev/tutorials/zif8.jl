#
# Images of the ZIF-8 tutorial (docs/src/tutorials/zif8.md): the initial placement (the
# `check_initial_placement` line of the input uncommented) and the packed system, seen
# along y, a slice of one unit cell (17 Å) in y, with whole molecules. Run from the root of
# the repository, with the packmol executable built (juliac/build.jl):
#
#     julia --project=docs dev/tutorials/zif8.jl
#
include("render.jl")
using PDBTools

const exe = joinpath(@__DIR__, "..", "..", "juliac", "build", "bin", "packmol")
const files = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files", "zif8")
const images = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "images")
const work = mktempdir()

for f in readdir(files)
    cp(joinpath(files, f), joinpath(work, f); force = true)
end
input = read(joinpath(work, "zif8.inp"), String)
write(joinpath(work, "check_zif8.inp"), replace(input,
    "# check_initial_placement" => "check_initial_placement",
    "output zif8_ethanol_water.pdb" => "output check_zif8_ethanol_water.pdb",
))
for inp in ("check_zif8.inp", "zif8.inp")
    run(pipeline(`$exe -i $(joinpath(work, inp))`; stdout = devnull))
end

# The molecules (whole) with the center in the slice 0 ≤ y < 17.033, and the framework
# atoms in it, written with the cell of the slice
function slice(file)
    atoms = read_pdb(file)
    keep = Atom[]
    for mol in eachresidue(atoms)
        y = sum(at.y for at in mol) / length(mol)
        if resname(mol) == "ZIF"
            append!(keep, [at for at in mol if 0 <= at.y < 17.033])
        elseif 0 <= y < 17.033
            append!(keep, collect(mol))
        end
    end
    s = Structure(keep)
    s.unitcell = UnitCell(51.099, 17.033, 80.0, 90.0, 90.0, 90.0)
    out = tempname() * ".pdb"
    write_pdb(out, s)
    return out
end

reps = [
    (selection = "resname ZIF and not element H", style = "Licorice 0.25 12 12", color = "ColorID 2"),
    (selection = "resname ZIF and element Zn", style = "VDW 0.6 20", color = "ColorID 6"),
    (selection = "resname ETOH", style = "VDW 1.0 12", color = "ColorID 3"),
    (selection = "resname TIP3", style = "VDW 1.0 12", color = "ColorID 0"),
]
for (file, image) in (("check_zif8_ethanol_water.pdb", "zif8_check.webp"), ("zif8_ethanol_water.pdb", "zif8.webp"))
    render(slice(joinpath(work, file)), joinpath(images, image);
        reps, rotate = [(:x, -90)], cell = :unitcell, size = (1400, 1800), width = 800)
end
