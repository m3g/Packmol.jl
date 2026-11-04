# Liquid crystals

Liquid crystals are made of elongated molecules, which keep a common orientation (the
*director*) while flowing like a liquid. Here, the molecule is 5CB
(4-cyano-4'-pentylbiphenyl), the classic room-temperature nematic liquid crystal,
packed in two arrangements: a nematic phase, with the molecules aligned along `z`, and a
twisted nematic cell, the arrangement of a pixel of a liquid crystal display, with the
director turning in the plane of the cell from one layer to the next. The orientations are
set by the [rotation constraint](@ref "Rotation constraint"), `constrain_rotation`, which
keeps the angle between an axis and a vector of the molecule within a range.

Files: [`5cb.pdb`](files/liquid_crystal/5cb.pdb), and the input files
[`nematic.inp`](files/liquid_crystal/nematic.inp) and
[`twisted.inp`](files/liquid_crystal/twisted.inp).

## The molecule

`5cb.pdb` is 5CB with its pentyl chain in the extended (all-trans) conformation (built
from its SMILES string, `CCCCCc1ccc(cc1)-c1ccc(cc1)C#N`, with Open Babel, and minimized
with the MMFF94 force field). Atom 1 is the end of the alkyl chain, and atom 19 the
nitrogen of the cyano group: the vector from atom 1 to atom 19, 15.8 Å long, is the long
axis of the molecule, and it is the vector given to `constrain_rotation`, in its `atoms`
(from) and `to` lines.

The molecules of a nematic phase are aligned along the director, but as many point one
way as the other: the phase is apolar. So each arrangement has two structure blocks, one
with the vector of the molecules along the director, and one with it opposite.

## The numbers of molecules

The density of 5CB is 1.02 g/mL near room temperature, and its molar mass 249.35 g/mol:
0.00246 molecules per Å³. The nematic phase fills a periodic cube of 50 Å
(125000 Å³) with 308 molecules, 154 in each direction. Each layer of the twisted cell is
50 × 50 × 15 Å (37500 Å³), with 92 molecules, 46 in each direction.

## Nematic phase

### Command line

```@eval
Main.input_file("files/liquid_crystal/nematic.inp")
```

`constrain_rotation z 0. 20.` keeps the vector from atom 1 to atom 19 within 20° of
`+z`, and `constrain_rotation z 180. 20.` within 20° of `-z`. Run it with:

```
packmol nematic.inp
```

It packs in a few seconds. The cone of 20° sets the order of the phase: the nematic
order parameter, ``S = \langle 3\cos^2\theta - 1 \rangle / 2`` (``\theta`` the angle
of a molecule with the director), is 0.89 in the packed system, while 5CB has about 0.6
near room temperature: a cone of 45° gives 0.59.

### Julia

```@setup liquid_crystal
cp(joinpath("files", "liquid_crystal", "5cb.pdb"), "5cb.pdb"; force = true)
```

```@example liquid_crystal
using Packmol, PDBTools

# The vector from atom 1 (end of the chain) to atom 19 (cyano N) within 20° of +z and -z
up = structure_type("5cb.pdb"; number = 154,
    constrain_rotation = Dict(:z => (0.0, 20.0, [1], [19])),
)
down = structure_type("5cb.pdb"; number = 154,
    constrain_rotation = Dict(:z => (180.0, 20.0, [1], [19])),
)
nematic = PackmolSystem([up, down];
    output = "nematic.pdb",
    tolerance = 2.0,
    unitcell = [50. 0. 0.; 0. 50. 0.; 0. 0. 50.], # centered at the origin
    # check_initial_placement = true,
)
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(nematic)
end # hide
nematic.status
```

The packed system, with the `visualize` function of
[PDBTools.jl](https://github.com/m3g/PDBTools.jl) (drag to rotate, scroll to zoom), the
molecules pointing up in red and those pointing down in blue:

```@example liquid_crystal
visualize(get_atoms(nematic),
    "chain A" => (style = :sticks, color = "red"),
    "chain B" => (style = :sticks, color = "blue");
    unitcell = UnitCell(nematic.unitcell),
    unitcell_origin = -sum(eachcol(nematic.unitcell)) / 2,
)
```

## Twisted nematic cell

In a twisted nematic cell, the molecules lie in the plane of the cell, and the director
turns by 90° from one wall of the cell to the other: from `x` at the bottom to `y` at the
top. Here, the cell has four layers of 15 Å, with the director at 0°, 30°, 60° and 90° from
`x`. The cell is periodic in `x` and `y`; along `z`, 10 Å of empty space above the top
layer separate it from the bottom one, in the periodic image of the cell.

The direction of a molecule in a plane takes more than one angle: at 30° from `x`, the
vector is also at 60° from `y` (at -30° from `x`, it would be at 120°), and at 90° from
`z`. So each structure block has three rotation constraints, one for each axis: for a
director at an angle φ from `x`, the vector is at φ from `x`, 90° - φ from `y`, and 90°
from `z`, and the opposite one at 180° - φ, 90° + φ and 90°.

### Command line

```@raw html
<details><summary>The input file <code>twisted.inp</code> (click to expand)</summary>
```
```@eval
Main.input_file("files/liquid_crystal/twisted.inp")
```
```@raw html
</details><br>
```

```
packmol twisted.inp
```

It packs in about 15 seconds with 14 threads.

### Julia

The same input, with the structure types of the layers built in a loop:

```@example liquid_crystal
# The vector from atom 1 to atom 19 at angles (x, y, z) from the axes, within 15° of
# the angles with x and y, and 10° of the plane
function layer(z1, z2, (ax, ay, az))
    return structure_type("5cb.pdb"; number = 46,
        constraints = [AbovePlane([0., 0., 1.], z1), BelowPlane([0., 0., 1.], z2)],
        constrain_rotation = Dict(
            :x => (ax, 15.0, [1], [19]),
            :y => (ay, 15.0, [1], [19]),
            :z => (az, 10.0, [1], [19]),
        ),
    )
end
structure_types = Packmol.StructureType{3,Float64}[]
for (k, φ) in enumerate((0.0, 30.0, 60.0, 90.0))
    z1 = -35.0 + 15 * (k - 1)
    push!(structure_types, layer(z1, z1 + 15, (φ, 90 - φ, 90.0)))
    push!(structure_types, layer(z1, z1 + 15, (180 - φ, 90 + φ, 90.0)))
end
twisted = PackmolSystem(structure_types;
    output = "twisted.pdb",
    tolerance = 2.0,
    unitcell = [50. 0. 0.; 0. 50. 0.; 0. 0. 70.], # centered at the origin
)
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(twisted)
end # hide
twisted.status
```

The layers, from the bottom to the top, in red, orange, green and blue (each structure
type is written with a chain identifier of its own, from `A` to `H`): seen from above,
the director turns from one layer to the next.

```@example liquid_crystal
visualize(get_atoms(twisted),
    "chain A or chain B" => (style = :sticks, color = "red"),
    "chain C or chain D" => (style = :sticks, color = "orange"),
    "chain E or chain F" => (style = :sticks, color = "green"),
    "chain G or chain H" => (style = :sticks, color = "blue");
    unitcell = UnitCell(twisted.unitcell),
    unitcell_origin = -sum(eachcol(twisted.unitcell)) / 2,
)
```
