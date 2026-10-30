# Gases in a zeolite

Methane and carbon dioxide adsorbed in the cages of a zeolite: chabazite (CHA), whose
cages, connected by windows of eight oxygen atoms (3.8 Å), are those of the zeolites used
to separate carbon dioxide from methane. The zeolite is a fixed molecule filling its
periodic cell, which is triclinic, and the gas molecules are packed in its free space.

An interactive view of the packed system is at the end of the [Julia](#Julia) section.

Files: [`cha.pdb`](files/zeolite/cha.pdb), [`methane.pdb`](files/zeolite/methane.pdb),
[`co2.pdb`](files/zeolite/co2.pdb), and the input file [`zeolite.inp`](files/zeolite/zeolite.inp).

## The structure of the zeolite

The [IZA Database of Zeolite Structures](https://europe.iza-structure.org/IZA-SC/) provides
the structures of the zeolite frameworks as CIF files, with the atoms of the asymmetric
unit and the symmetry of the crystal. The full unit cell is obtained with
[Open Babel](https://openbabel.org):

```
obabel CHA.cif -O CHA_uc.pdb --fillUC strict
```

The unit cell of chabazite, in its hexagonal setting, has `a = b = 13.675` Å,
`c = 14.767` Å, and an angle of 120° between `a` and `b`, with 36 silicon and 72 oxygen
atoms, and 3 cages. [`cha.pdb`](files/zeolite/cha.pdb) has 2 × 2 × 2 unit cells (864
atoms, 24 cages), centered at the origin, built with the Julia code at the end of this
page.

## Command line

The input file:

```
#
# Methane and carbon dioxide in the cages of the chabazite zeolite (CHA):
# 2 x 2 x 2 hexagonal unit cells of the IZA database structure (24 cages), with
# one molecule of each gas per cage, in the triclinic periodic cell of the
# zeolite (a = b = 27.35 Å, c = 29.534 Å, gamma = 120 degrees).
#
tolerance 2.0
output zeolite.pdb
unitcell 27.350 27.350 29.534 90. 90. 120.

# The zeolite, kept at the coordinates of the file (centered at the origin, as
# the cell). The radius of the atoms (1.5 Å, here and in the gases) keeps the
# gases 3 Å away from the framework.
structure cha.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  radius 1.5
end structure

structure methane.pdb
  number 24
  radius 1.5
end structure

structure co2.pdb
  number 24
  radius 1.5
end structure
```

The periodic cell, given by `unitcell a b c alpha beta gamma` (see [Triclinic
box](../periodic_boundary_conditions/triclinic_box.md)), is the cell of the zeolite, so the
framework is continuous across the periodic boundaries. The cell of `unitcell` is centered
at the origin, as the coordinates of [`cha.pdb`](files/zeolite/cha.pdb), so that all the
atoms of the fixed zeolite are inside it. No constraint is needed for the gases: the cell
confines them, and the distances to the atoms of the fixed framework keep them in the
cages.

With the default `tolerance` of 2 Å, the atoms of the gases could get 2 Å from the
oxygen atoms of the framework, close enough to fit in rings too small for any molecule.
The `radius` of the atoms sets the minimum distance between atoms of different molecules
to the sum of their radii: 3 Å here, a typical contact distance between these atoms. The
`radius` keyword applies to the atoms of one structure, so the radii can be different for
each molecule (see [Structure blocks](../input_files/structure_blocks.md)).

### Check the initial approximation first

Add the `check` keyword to the input file, and run it as below. Packmol builds only the
initial approximation, in a fraction of a second, and writes it to the output file
(`zeolite.pdb`): look at it to see that the zeolite fills the periodic cell. When it does,
remove `check`, and run the packing.

### Packing

```
packmol zeolite.inp
```

The packing takes a fraction of a second. The output file has the `CRYST1` record of the
triclinic cell, ready for a simulation.

## Julia

```@setup zeolite
for file in ("cha.pdb", "methane.pdb", "co2.pdb")
    cp(joinpath("files", "zeolite", file), file; force = true)
end
```

```@example zeolite
using Packmol, PDBTools

cha = structure_type("cha.pdb";
    number = 1,
    fixed = ([0., 0., 0.], [0., 0., 0.]),
    radius = 1.5,
)
methane = structure_type("methane.pdb"; number = 24, radius = 1.5)
co2 = structure_type("co2.pdb"; number = 24, radius = 1.5)
system = PackmolSystem([cha, methane, co2];
    output = "zeolite.pdb",
    tolerance = 2.0,
    unitcell = UnitCell(27.350, 27.350, 29.534, 90., 90., 120.), # centered at the origin
)
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
system.status
```

`UnitCell` (from [PDBTools.jl](https://github.com/m3g/PDBTools.jl)) builds the matrix of the
cell vectors from its lattice parameters. As with the input file, check the initial
approximation first, with `PackmolSystem(...; check = true)`.

The packed system, with the `visualize` function of PDBTools.jl, and the edges of the
triclinic cell (whose corner is at `-(a + b + c)/2`, since it is centered at the origin):
the framework of the zeolite, methane and carbon dioxide (drag to rotate, scroll to zoom):

```@example zeolite
visualize(get_atoms(system);
    unitcell = UnitCell(system.unitcell),
    unitcell_origin = -sum(eachcol(system.unitcell)) / 2,
)
```

### Building the supercell

The 2 × 2 × 2 cells of the zeolite, from the unit cell given by Open Babel, centered at
the origin:

```julia
using PDBTools

uc = read_pdb("CHA_uc.pdb")
cell = Matrix{Float64}(UnitCell(13.675, 13.675, 14.767, 90, 90, 120))
A, B, C = cell[:, 1], cell[:, 2], cell[:, 3] # the cell vectors
center = A + B + C # the center of the supercell
atoms = Atom{Nothing}[]
for i in 0:1, j in 0:1, k in 0:1, at in uc
    el = element(at) # Si or O
    x = position(at) + i * A + j * B + k * C - center
    push!(atoms, Atom(; name = el, resname = "CHA", resnum = 1, chain = "A",
        x = x[1], y = x[2], z = x[3], pdb_element = el,
    ))
end
supercell = Structure(atoms)
supercell.unitcell = UnitCell(2 * 13.675, 2 * 13.675, 2 * 14.767, 90, 90, 120)
write_pdb("cha.pdb", supercell)
```

The gases are small enough to be written by hand: methane, with C–H bonds of 1.09 Å in a
tetrahedron, and carbon dioxide, with C=O bonds of 1.16 Å (see the files).
