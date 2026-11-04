# Liquid crystals

Liquid crystals are made of elongated molecules, which keep a common orientation (the
*director*) while flowing like a liquid. Here, the molecule is 5CB
(4-cyano-4'-pentylbiphenyl), the classic room-temperature nematic liquid crystal,
packed in two arrangements: a nematic phase, with the molecules aligned along `z`, and a
cholesteric (twisted) phase, with the director turning in the xy plane from one layer to
the next, along a helix. The orientations are
set by the [rotation constraint](@ref "Rotation constraint"), `constrain_rotation`, which
keeps the angle between an axis and a vector of the molecule within a range.

Files: [`5cb.pdb`](files/liquid_crystal/5cb.pdb), and the input files
[`nematic.inp`](files/liquid_crystal/nematic.inp) and
[`cholesteric.inp`](files/liquid_crystal/cholesteric.inp).

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
(125000 Å³) with 308 molecules, 154 in each direction. Each layer of the cholesteric
phase is 50 × 50 × 15 Å (37500 Å³), with 92 molecules, 46 in each direction.

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

## Cholesteric phase

In a cholesteric (chiral nematic) phase, the director turns along a helix: the molecules
of each layer lie in the xy plane, aligned along a director that turns from one layer to
the next. 5CB forms it when a chiral molecule is added to it; in a twisted nematic cell,
the arrangement of a pixel of a liquid crystal display, the twist is instead imposed by the
walls of the cell. Here, the cell has four layers of 15 Å, with the director at 0°, 45°,
90° and 135° from `x`: the next layer, through the periodic boundary along `z`, would be
at 180°, the same as 0° for an apolar director, so the twist continues through the
boundary, along a helix of a pitch of 120 Å (a full turn of the director).

The direction of a molecule in a plane takes more than one angle: at 45° from `x`, the
vector is also at 45° from `y` (at -45° from `x`, it would be at 135° from `y`), and at
90° from `z`. So each structure block has three rotation constraints, one for each axis:
for a vector at an angle φ from `x` in the plane, `(cos φ, sin φ, 0)`, its angles are
`acos(cos φ)` with `x`, `acos(sin φ)` with `y`, and 90° with `z`. The molecules pointing
the other way are at φ + 180°.

### Command line

```@raw html
<details><summary>The input file <code>cholesteric.inp</code> (click to expand)</summary>
```
```@eval
Main.input_file("files/liquid_crystal/cholesteric.inp")
```
```@raw html
</details><br>
```

```
packmol cholesteric.inp
```

It packs in about 10 seconds with 14 threads.

### Julia

The same input, with the structure types of the layers built in a loop:

```@example liquid_crystal
# The vector from atom 1 to atom 19 at the angle φ from x, in the xy plane: within 15°
# of its angles with x and y, and 10° of the plane
function layer(z1, z2, φ)
    return structure_type("5cb.pdb"; number = 46,
        constraints = [AbovePlane([0., 0., 1.], z1), BelowPlane([0., 0., 1.], z2)],
        constrain_rotation = Dict(
            :x => (acosd(cosd(φ)), 15.0, [1], [19]),
            :y => (acosd(sind(φ)), 15.0, [1], [19]),
            :z => (90.0, 10.0, [1], [19]),
        ),
    )
end
structure_types = Packmol.StructureType{3,Float64}[]
for (k, φ) in enumerate((0.0, 45.0, 90.0, 135.0))
    z1 = -30.0 + 15 * (k - 1)
    push!(structure_types, layer(z1, z1 + 15, φ))
    push!(structure_types, layer(z1, z1 + 15, φ + 180)) # antiparallel
end
cholesteric = PackmolSystem(structure_types;
    output = "cholesteric.pdb",
    tolerance = 2.0,
    unitcell = [50. 0. 0.; 0. 50. 0.; 0. 0. 60.], # centered at the origin
)
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(cholesteric)
end # hide
cholesteric.status
```

The layers, from the bottom to the top, in red, orange, green and blue (each structure
type is written with a chain identifier of its own, from `A` to `H`): seen from above,
the director turns by 45° from one layer to the next.

```@example liquid_crystal
visualize(get_atoms(cholesteric),
    "chain A or chain B" => (style = :sticks, color = "red"),
    "chain C or chain D" => (style = :sticks, color = "orange"),
    "chain E or chain F" => (style = :sticks, color = "green"),
    "chain G or chain H" => (style = :sticks, color = "blue");
    unitcell = UnitCell(cholesteric.unitcell),
    unitcell_origin = -sum(eachcol(cholesteric.unitcell)) / 2,
)
```
