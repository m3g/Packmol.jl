# Surfaces given by points

```
within surface     file distance [offset]
not within surface file distance [offset]
inside surface     file [offset]
outside surface    file [offset]
```

A surface, or a curve, given by points in a file: the atoms are kept within `distance`
(Å) of it, or farther than `distance` from it (`not within`). It is how membranes and
compartments of any shape are built: the points of a triangulated surface, of a
[TS2CG](https://github.com/marrink-lab/TS2CG) model, of a shape drawn in a 3D editor or
sampled from any function.

The file is relative to the directory of the input file, and holds a point per line,
`x y z` or `x y z nx ny nz` (Å, with the normal of the surface at the point; lines
starting with `#` are comments), or it is a point file of TS2CG (the `.dat` files of its
`PLM` module, as `point/OuterBM.dat`, in nm: they are converted to Å).

**With normals**, the distance `h` to the surface is signed: positive on the side to
which the normals point. Then `offset` shifts the surface along its normals, and the
constraint is `|h - offset| ≤ distance` (`within`) or `|h - offset| ≥ distance`
(`not within`). A single file of points of the middle of a bilayer, for example, defines
where the heads of the lipids of both leaflets go (offsets of ±19 Å, for POPC), where the
ends of their chains go (offset 0), and where the solvent goes (not within 24 Å of it).
**Without normals**, `h` is the distance to the nearest point, and `offset` must be zero:
for curves, or when the side of the surface does not matter.

The signed distance is smooth: near the surface it is the average of the distances to
the tangent planes of the points within a smoothing radius (an implicit moving least
squares surface), so that the points need not be dense, as long as they sample the
shape: the default radius is 2.5 times the typical distance between neighboring points.
On curved surfaces, the smoothing makes the surface slightly flatter: about `r²/(10R)`
inwards, for a smoothing radius `r` and a radius of curvature `R` (0.2 Å, for `r` = 14 Å
and `R` = 100 Å).

`not within` is the cheaper constraint: an atom with no point of the surface within the
smoothing radius satisfies it (if `distance + |offset|` is not larger than the radius),
without further search.

**Inside and outside.** For a closed surface whose normals point outwards (as those of
the outer membrane of a vesicle of TS2CG, or of any closed shape), `inside surface` keeps
the atoms inside it, `h ≤ offset`, and `outside surface` outside it, `h ≥ offset`. With
`offset`, the region is bounded by the surface parallel to it at that distance along the
normals: `inside surface cell.dat -5.0` keeps the atoms 5 Å within the surface,
`outside surface cell.dat 3.0` 3 Å beyond it. These constraints need normals, and the
normals decide the side: they must point out of the enclosed region everywhere, and the
surface must have no holes wider than the smoothing radius. Far from the surface, the
side is that of a precomputed cell of the grid of the points, so that the solvent filling
a large compartment costs little. For the initial placement and the volume tools, the
region inside the surface is bounded by the box of its points.

In Julia code, the surface is a [`PointSurface`](@ref), built from vectors of points and
normals or read from a file with [`read_surface`](@ref), and the constraints are
[`WithinSurface`](@ref), [`NotWithinSurface`](@ref), [`InsideSurface`](@ref) and
[`OutsideSurface`](@ref). Several constraints can share a
surface; in input files, constraints naming the same file share it.

```@docs
PointSurface
read_surface
surface_distance
WithinSurface
NotWithinSurface
InsideSurface
OutsideSurface
```

## Example

A bilayer of a coarse-grained lipid on a torus given by points with normals, and some
water around it, not within 12 Å of the middle of the bilayer. The lipid is a chain of 5
beads, 8 Å long: the head (bead 1) goes at 9 Å from the middle surface, on either side,
and the end of the chain (bead 5) near it. The numbers of lipids follow from the areas of
the two leaflets (about 60 Å² per lipid):

```
structure lipid.pdb
  number 530
  atoms 1
    within surface torus.dat 1.5 9.0
  end atoms
  atoms 5
    within surface torus.dat 1.5
  end atoms
end structure

structure lipid.pdb
  number 115
  atoms 1
    within surface torus.dat 1.5 -9.0
  end atoms
  atoms 5
    within surface torus.dat 1.5
  end atoms
end structure

structure water.pdb
  number 1000
  not within surface torus.dat 12.0
  inside box -60. -60. -25. 60. 60. 25.
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl). The points
of the torus (radii of 35 and 14 Å around the `z` axis) are generated here, with their
normals:

```@setup surface
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
# A 5-bead lipid along z: the head at the top
open("lipid.pdb", "w") do io
    for i in 1:5
        println(io, "ATOM  $(lpad(i, 5))  C$i  LIP     1       0.000   0.000$(lpad(string(round(10.0 - 2.0 * i; digits = 3)), 8))  1.00  0.00")
    end
end
```

```@example surface
using Packmol, PDBTools
# Points of the torus, about 1.5 Å apart, and the normals
R = 35.0 # the radius of the circle of the centers of the tube
r = 14.0 # the radius of the tube
points = NTuple{3,Float64}[]
normals = NTuple{3,Float64}[]
for u in range(0, 2π; length = 200)[1:end-1], v in range(0, 2π; length = 70)[1:end-1]
    n = (cos(v) * cos(u), cos(v) * sin(u), sin(v))
    push!(points, (R * cos(u) + r * n[1], R * sin(u) + r * n[2], r * n[3]))
    push!(normals, n)
end
torus = PointSurface(points, normals)
outer = structure_type("lipid.pdb"; number = 530,
    atom_constraints = [[1] => WithinSurface(torus, 1.5; offset = 9.0), [5] => WithinSurface(torus, 1.5)],
)
inner = structure_type("lipid.pdb"; number = 115,
    atom_constraints = [[1] => WithinSurface(torus, 1.5; offset = -9.0), [5] => WithinSurface(torus, 1.5)],
)
water = structure_type("water.pdb"; number = 1000,
    constraints = [NotWithinSurface(torus, 12.0), InsideBox([-60., -60., -25.], [60., 60., 25.])],
)
system = PackmolSystem([outer, inner, water]; tolerance = 2.0, output = "torus.pdb")
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system);
    rotate = [(:z, 30), (:x, -70)],
)
```

## Example: inside and outside

Water inside the same torus, and outside it, at least 3 Å away from its surface, in a
box:

```
structure water.pdb
  number 400
  inside surface torus.dat
end structure

structure water.pdb
  number 1500
  outside surface torus.dat 3.0
  inside box -60. -60. -25. 60. 60. 25.
end structure
```

```@example surface
inside = structure_type("water.pdb"; number = 400, constraints = [InsideSurface(torus)])
outside = structure_type("water.pdb"; number = 1500,
    constraints = [OutsideSurface(torus; offset = 3.0), InsideBox([-60., -60., -25.], [60., 60., 25.])],
)
system = PackmolSystem([inside, outside]; tolerance = 2.0, output = "torus_water.pdb")
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system);
    rotate = [(:z, 30), (:x, -70)],
)
```
