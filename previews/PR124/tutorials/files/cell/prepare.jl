#
# Packmol inputs for the minimal cell JCVI-Syn3A from the published model of Bentopy
# (Westendorp et al., Protein Science 2026, doi:10.1002/pro.70480; data: Zenodo record
# 17524049, doi:10.5281/zenodo.17524049, CC-BY 4.0).
#
#   julia prepare.jl path/to/zenodo/cell
#
# `path/to/zenodo/cell` is the `cell` directory of the Zenodo archive (`zenodo.tar.gz`).
# Writes, in the current directory:
#   structures/     the molecules of the cell (PDB, Å), from Bentopy's structures
#   membrane.cif    the membrane (with its proteins), from Bentopy's model, fixed
#   chromosome.cif  the chromosome, from Bentopy's model, fixed
#   cell.inp        the proteins, RNAs and metabolites of Bentopy's input (input.bent),
#                   in the cytosol, placed by voxel_place and packed
#   w.pdb, gly.pdb  the solvent: a Martini water bead, and a two-bead placeholder of
#                   glycerol (for solvation.inp)
#
using PDBToolsBase: Atom, read_pdb, write_pdb, write_mmcif
using LinearAlgebra: norm
using Printf: @sprintf
using Statistics: mean

const bentopy = ARGS[1]

# A .gro file as atoms (nm to Å)
function read_gro(file)
    io = open(file)
    readline(io)
    n = parse(Int, strip(readline(io)))
    atoms = Vector{Atom}(undef, n)
    for i in 1:n
        l = readline(io)
        atoms[i] = Atom(; index = i, index_pdb = i, resnum = parse(Int, l[1:5]),
            resname = strip(l[6:10]), name = strip(l[11:15]),
            x = 10parse(Float64, l[21:28]), y = 10parse(Float64, l[29:36]), z = 10parse(Float64, l[37:44]))
    end
    close(io)
    return atoms
end

# The membrane and the chromosome, fixed; the center of the cell: the center of the lipids
membrane = read_gro(joinpath(bentopy, "membrane.gro"))
chromosome = read_gro(joinpath(bentopy, "chromosome.gro"))
lipids = ("CHOL", "CDL2", "DPSM", "POPC", "POPE", "POPG", "DPPC", "POPS", "DOPC", "DPPE", "PPCS")
center = mean([Float64(a.x), Float64(a.y), Float64(a.z)] for a in membrane if a.resname in lipids)   # Float64: 7.8M atoms
write_mmcif("membrane.cif", membrane)
write_mmcif("chromosome.cif", chromosome)

# The cytosol: within 950 Å of the center (the reach of the cytosol of Bentopy's model; the
# membrane bounds it). The peripheral proteins (Bentopy: within 20 nm of the membrane): the
# atom closest to their center beyond 630 Å.
R_cytosol = 950.0
R_peripheral = 630.0
sphere(R) = @sprintf("sphere %.4f %.4f %.4f %.1f", center..., R)

mkpath("structures")
open("cell.inp", "w") do io
    println(io, """
    #
    # The minimal cell JCVI-Syn3A, from the model of Bentopy (Westendorp et al., Protein
    # Science 2026): its membrane and chromosome fixed, and its proteins, RNAs and
    # metabolites placed by the voxel placement of Bentopy (voxel_place) and packed.
    #
    # check_initial_placement
    tolerance 4.0
    seed 1
    output cell.cif

    # The voxel placement of Bentopy, with its max-tries-mult for this model
    voxel_place 100000

    structure membrane.cif
      number 1
      fixed 0. 0. 0. 0. 0. 0.
    end structure

    structure chromosome.cif
      number 1
      fixed 0. 0. 0. 0. 0. 0.
    end structure
    """)
    for l in eachline(joinpath(bentopy, "input.bent"))
        m = match(r"^(\S+)\s+(\d+)\s+from\s+\"([^\"]+)\"\s+in\s+(\S+)", l)
        isnothing(m) && continue
        name, number, path, compartment = m.captures
        parse(Int, number) > 0 || continue
        name = replace(split(name, ':')[1], r"[^A-Za-z0-9_-]" => "_")
        file = "structures/$name.pdb"
        atoms = endswith(path, ".gro") ? read_gro(joinpath(bentopy, path)) : read_pdb(joinpath(bentopy, path))
        isfile(file) || write_pdb(file, atoms)
        println(io, "structure $file\n  number $number\n  inside $(sphere(R_cytosol))")
        if compartment == "peripheral"
            x = [[a.x, a.y, a.z] for a in atoms]
            icenter = argmin(norm.(x .- Ref(mean(x))))
            println(io, "  atoms $icenter\n    outside $(sphere(R_peripheral))\n  end atoms")
        end
        println(io, "end structure")
    end
end

# The solvent
write("w.pdb", "HETATM    1  W   W       1       0.000   0.000   0.000  1.00  0.00           W\nEND\n")
write("gly.pdb", """
    HETATM    1  G1  GLY     1       0.000   0.000   0.000  1.00  0.00           C
    HETATM    2  G2  GLY     1       4.000   0.000   0.000  1.00  0.00           C
    END
    """)
println("center of the cell: ", center)
