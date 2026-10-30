#
# Single-molecule structure files for the tutorials, extracted from the systems of
# ComplexMixturesExamples (https://github.com/m3g/ComplexMixturesExamples), whose
# equilibrated structures are in ~/Documents/ComplexMixturesExamples/scripts.
#
using PDBTools
const cme = expanduser("~/Documents/ComplexMixturesExamples/scripts")
const files = joinpath(@__DIR__, "..", "..", "docs", "src", "tutorials", "files")

# The first molecule of `resname` (a residue), or the atoms selected by `sel`
function extract(src, dir, name, sel)
    atoms = read_pdb(joinpath(cme, src), sel)
    if any(at -> resname(at) in ("TIP3", "GLYC", "DMF", "ETOH", "POPC"), atoms)
        first_residue = resnum(atoms[1])
        seg = segname(atoms[1])
        atoms = filter(at -> resnum(at) == first_residue && segname(at) == seg, atoms)
    end
    # Centered at the origin, in a chain of its own
    c = sum(position(at) for at in atoms) / length(atoms)
    for at in atoms
        at.x, at.y, at.z = position(at) - c
        at.chain = "A"
        at.segname = ""
    end
    mkpath(joinpath(files, dir))
    write_pdb(joinpath(files, dir, name), atoms)
    println(rpad(joinpath(dir, name), 40), length(atoms), " atoms")
end

extract("example1/system.pdb", "protein_glycerol", "protein.pdb", "protein")
extract("example1/system.pdb", "protein_glycerol", "glycerol.pdb", "resname GLYC")
extract("example1/system.pdb", "protein_glycerol", "water.pdb", "resname TIP3")
extract("example2/equilibrated.pdb", "polymer_dmf", "polyacrylamide.pdb", "resname ACR FACR LACR")
extract("example2/equilibrated.pdb", "polymer_dmf", "dmf.pdb", "resname DMF")
extract("example3/equilibrated.pdb", "membrane", "popc.pdb", "resname POPC")
extract("example3/equilibrated.pdb", "membrane", "ethanol.pdb", "resname ETOH")
extract("example3/equilibrated.pdb", "membrane", "water.pdb", "resname TIP3")
