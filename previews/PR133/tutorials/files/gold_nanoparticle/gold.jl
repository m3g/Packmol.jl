#
# The gold core of the nanoparticle (gold.pdb): the atoms of the fcc lattice of gold
# (a = 4.078 Å) within 15 Å of an atom, centered at the origin.
#
#   julia gold.jl gold.pdb
#
using PDBTools

a = 4.078
R = 15.0
basis = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]
n = ceil(Int, R / a) + 1
atoms = Atom{Nothing}[]
for i in -n:n, j in -n:n, k in -n:n, (bx, by, bz) in basis
    x, y, z = a .* (i + bx, j + by, k + bz)
    x^2 + y^2 + z^2 <= R^2 || continue
    push!(atoms, Atom(; name = "AU", resname = "AU", resnum = 1, chain = "A",
        x, y, z, pdb_element = "Au"))
end
write_pdb(ARGS[1], atoms)
println(length(atoms), " gold atoms, radius ", round(maximum(sqrt(at.x^2 + at.y^2 + at.z^2) for at in atoms); digits = 2), " Å")
