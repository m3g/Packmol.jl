# Water and urea

A periodic box of 1000 water molecules and 400 urea molecules: the simplest kind of
system, a homogeneous liquid, in which every molecule can be anywhere in the cell.

```@setup mixture_view
using Packmol, PDBTools
system = Main.pack_tutorial("mixture", "mixture.inp")
```

```@example mixture_view
visualize(get_atoms(system), # hide
    "water" => (style = :sticks, color = :element, water = true), # hide
    "not water" => (style = :ballandstick, color = :element); # hide
    unitcell = UnitCell(system.unitcell), # hide
    unitcell_origin = system.unitcell_center - sum(eachcol(system.unitcell)) / 2, # hide
    rotate = [(:z, 30), (:x, -70)], # hide
) # hide
```

*Water (sticks) and urea (balls and sticks), with the edges of the periodic box (drag to rotate, scroll to zoom).*

Files: [`water.pdb`](files/mixture/water.pdb), [`urea.pdb`](files/mixture/urea.pdb) and
the input file [`mixture.inp`](files/mixture/mixture.inp).

## The periodic box

The input file:

```@eval
Main.input_file("files/mixture/mixture.inp")
```

`pbc 40. 40. 40.` sets a periodic box from the origin to `(40, 40, 40)` (see [Orthorhombic
box](../periodic_boundary_conditions/orthorhombic_box.md)). The box confines the
molecules, so they need no constraint, and the distances between them are computed
across its faces: the molecules near a face are kept `tolerance` away from those near the
opposite face, and the box can be used in a simulation with periodic boundary conditions
as it is, with no gaps or overlaps at the faces. The output file has the `CRYST1` record
of the box.

Run it with:

```
packmol mixture.inp
```

The packed system is written to `mixture.pdb`. A system
this simple needs no check of the initial placement, but the input has the
`check_initial_placement` line, commented, as all the tutorials (see
[Checking the initial placement](@ref check_initial_placement)).

## Julia script

Each `structure ... end structure` block is a [`structure_type`](@ref), and the global
keywords (`tolerance`, `output`, `pbc`) are keyword arguments of [`PackmolSystem`](@ref):

```julia
using Packmol

water = structure_type("water.pdb"; number = 1000)
urea = structure_type("urea.pdb"; number = 400)

system = PackmolSystem([water, urea];
    output = "mixture.pdb",
    tolerance = 2.0,
    unitcell = [40. 0. 0.; 0. 40. 0.; 0. 0. 40.], # centered at the origin
    # check_initial_placement = true, # uncomment to only place the molecules
)
packmol(system)
```

`unitcell` is the matrix of the cell vectors (as columns); the cell is centered at the
origin, unless `unitcell_center` is given (`unitcell_center = [20., 20., 20.]` is the box
of `pbc 40. 40. 40.`). The output format follows the extension of the output file, so
`filetype` is not needed.

The packed coordinates are also available without reading the output file:
`get_atoms(system)` returns them as a vector of atoms.
