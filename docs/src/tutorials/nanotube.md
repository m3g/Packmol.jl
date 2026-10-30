# Carbon nanotube filled with water

A carbon nanotube filled with water, in water. The nanotube is a fixed molecule along the
`z` axis, crossing the periodic boundary of the box, so that it is infinite; the water
fills the inside of the tube, and the box around it.

Files: [`cnt.pdb`](files/nanotube/cnt.pdb), [`water.pdb`](files/nanotube/water.pdb), and
the input file [`nanotube.inp`](files/nanotube/nanotube.inp).

The nanotube is an armchair (12,12) tube, with a radius of 8.14 Å, and 49.2 Å long (20
periods of 2.46 Å along its axis, 960 carbon atoms), built with the Julia code at the end
of this page. Its length is the length of the box along `z`, so that the tube is continuous
across the periodic boundary.

## The numbers of molecules

The atoms of water can get no closer than `tolerance` (2 Å) to the atoms of the nanotube,
so inside the tube they are within a radius of about 6 Å. A cylinder of this radius and
49.2 Å long has 5560 Å³, which would hold about 185 molecules at 1 g/mL; water confined
in a tube this narrow is less dense, and 150 molecules are packed. Outside the tube,
beyond a radius of 10.5 Å, the box of 40 × 40 × 49.2 Å has 61700 Å³, which holds 2060
molecules at 1 g/mL (0.0334 molecules per Å³).

## Command line

The input file:

```
#
# A carbon nanotube (armchair, (12,12), radius 8.14 Å) filled with water, in
# water (1 g/mL), in a periodic box. The nanotube is along z, and crosses the
# periodic boundary: it is infinite.
#
tolerance 2.0
output nanotube.pdb
pbc -20. -20. -24.6 20. 20. 24.6

# The nanotube, kept at the coordinates of the file
structure cnt.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
end structure

# Water inside the tube. The cylinders start at z = -100 with a length of 200,
# longer than the box: they cross the periodic boundary along z.
structure water.pdb
  number 150
  inside cylinder 0. 0. -100. 0. 0. 1. 6. 200.
end structure

# Water outside the tube
structure water.pdb
  number 2060
  outside cylinder 0. 0. -100. 0. 0. 1. 10.5 200.
end structure
```

The cylinders (`center` of one cap, `axis`, `radius` and `length`) are longer than the box:
with periodic boundary conditions, the atoms are evaluated in the box, where only the side
of the cylinders matters.

### Check the initial approximation first

Add the `check` keyword to the input file, and run it as below. Packmol builds only the
initial approximation, in a fraction of a second, and writes it to the output file
(`nanotube.pdb`): see it at the end of the [Julia](#Julia) section. When the water is in
the tube and around it, remove `check`, and run the packing.

### Packing

```
packmol nanotube.inp
```

The packing takes a fraction of a second.

## Julia

```@setup nanotube
for file in ("cnt.pdb", "water.pdb")
    cp(joinpath("files", "nanotube", file), file; force = true)
end
```

```@example nanotube
using Packmol, PDBTools

cnt = structure_type("cnt.pdb"; number = 1, fixed = ([0., 0., 0.], [0., 0., 0.]))
water_inside = structure_type("water.pdb";
    number = 150,
    constraints = [InsideCylinder([0., 0., -100.], [0., 0., 1.], 6., 200.)],
)
water_outside = structure_type("water.pdb";
    number = 2060,
    constraints = [OutsideCylinder([0., 0., -100.], [0., 0., 1.], 10.5, 200.)],
)
system = PackmolSystem([cnt, water_inside, water_outside];
    output = "nanotube.pdb",
    tolerance = 2.0,
    unitcell = [40. 0. 0.; 0. 40. 0.; 0. 0. 49.2], # centered at the origin
    check = true, # first, only the initial approximation
)
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
system.status
```

The initial approximation, with the `visualize` function of
[PDBTools.jl](https://github.com/m3g/PDBTools.jl) (drag to rotate, scroll to zoom): the
water is in the tube and around it, but molecules may overlap:

```@example nanotube
show_system(system) = visualize(get_atoms(system),
    "resname CNT" => (style = :sticks, color = "gray"),
    "water" => (style = :sticks, water = true);
    unitcell = UnitCell(system.unitcell),
    unitcell_origin = -sum(eachcol(system.unitcell)) / 2,
)
show_system(system)
```

Then the packing:

```@example nanotube
system.check = false
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
show_system(system)
```

### Building the nanotube

The nanotube is the rectangular cell of graphene (2.46 × 4.26 Å, with 4 carbon atoms) rolled
into a cylinder: for an armchair `(n, n)` tube, `n` cells along the circumference, with the
C–C bonds along it, and the tube periodic along its axis with a period of 2.46 Å:

```julia
using PDBTools

function build_nanotube(; n = 12, nz = 20, a = 2.46, d = 1.42)
    cell = [(0.0, 0.0), (a / 2, d / 2), (a / 2, 3d / 2), (0.0, 2d)]
    circumference = n * 3d
    radius = circumference / (2π)
    Lz = nz * a
    atoms = Atom{Nothing}[]
    for i in 0:nz-1, j in 0:n-1, (z, y) in cell
        θ = 2π * (j * 3d + y) / circumference
        push!(atoms, Atom(; name = "C", resname = "CNT", resnum = 1, chain = "A",
            x = radius * cos(θ), y = radius * sin(θ), z = i * a + z - Lz / 2,
            pdb_element = "C",
        ))
    end
    return atoms
end

write_pdb("cnt.pdb", build_nanotube())
```
