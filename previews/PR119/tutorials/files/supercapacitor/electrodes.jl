#
# The electrodes of the supercapacitor (electrodes.pdb): two blocks of graphite, of three
# graphene layers each (AB stacked, 3.35 Å apart), periodic in x and y (39.36 x 42.60 Å).
# The cathode has its layers at z = 1.0, 4.35 and 7.7 Å, the anode at z = 62.7, 66.05 and
# 69.4 Å.
#
#   julia electrodes.jl electrodes.pdb
#
using PDBTools

nx, ny = 16, 10                   # rectangular cells of graphene, of 2.46 x 4.26 Å
a, b = 2.46, 4.26
basis = [(0.0, 0.0), (1.23, 0.71), (1.23, 2.13), (0.0, 2.84)]   # C–C bonds of 1.42 Å
layers = [(1.0, 0.0), (4.35, 1.42), (7.70, 0.0),        # (z, shift in y of the B layers)
          (62.70, 0.0), (66.05, 1.42), (69.40, 0.0)]
atoms = Atom{Nothing}[]
for (il, (z, dy)) in enumerate(layers), i in 0:nx-1, j in 0:ny-1, (bx, by) in basis
    resname, chain = il <= 3 ? ("CAT", "A") : ("ANO", "B")    # cathode, anode
    push!(atoms, Atom(; name = "C", resname, resnum = 1, chain,
        x = i * a + bx, y = mod(j * b + by + dy, ny * b), z, pdb_element = "C"))
end
write_pdb(ARGS[1], atoms)
println(length(atoms), " atoms, cell of ", nx * a, " x ", ny * b, " Å")
