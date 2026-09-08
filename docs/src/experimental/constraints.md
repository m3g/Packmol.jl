# Constraints

!!! warning "Experimental"
    This documents the native Julia packing engine's constraint syntax,
    shared between [input files](input_files.md) and the
    [Julia API](julia_api.md). It tracks the original Fortran
    [Packmol](http://github.com/m3g/packmol) keywords where they exist, but
    is not yet a complete port — see [Home](index.md) for what's covered so
    far.

A constraint restricts where the atoms of a structure are allowed to be. In
an input file, each one is declared as a line inside a
`structure ... end structure` block, and takes the form:

```
<placement> <shape> <parameters...>
```

Multiple constraints can be given for the same structure — a molecule must
then satisfy all of them simultaneously (e.g. `inside box ...` together with
`outside sphere ...` carves a box with a spherical hole in it). By default a
constraint applies to every atom of the structure; in an input file, wrapping
a subset of lines in an `atoms <indices> ... end` block restricts a
constraint to just those atoms (see [Structure blocks](input_files.md) in
Input files).

## Periodic boundary conditions

```
pbc a b c
pbc xmin ymin zmin  xmax ymax zmax
pbc dodecahedral cx cy cz d
unitcell a b c alpha beta gamma
```

`pbc`/`unitcell` are not shape constraints declared inside a `structure`
block — they are global keywords that set up a periodic simulation cell (see
[Input files](input_files.md)). `pbc` gives an orthorhombic box, either as
side lengths (centered at the origin) or explicit min/max corners, or a
rhombic dodecahedron centered at `(cx,cy,cz)` (see below); `unitcell` gives
a general triclinic cell, CRYST1-style, also centered at the origin. Once
set, interatomic distances (and any explicit constraint below) are
evaluated using the periodic cell — atoms are wrapped to the image centered
on the cell before a constraint checks them.

### Dodecahedral box

`pbc dodecahedral cx cy cz d` sets up a rhombic dodecahedron cell of "size"
`d` centered at `(cx,cy,cz)` — internally just the triclinic cell `a = b =
c = d`, `α = β = 60°`, `γ = 90°` (the standard MD convention, matching
GROMACS's `editconf -bt dodecahedron -d d`), so it's written out (e.g. in
`CRYST1`) and behaves exactly like any other `unitcell`. For a given cell
volume, this shape gets a bigger minimum-image distance than a cube would
(`d` is that guaranteed minimum-image distance, in every direction) — about
70.7% of a cube's volume is enough to enclose the same sphere, which is why
it's the usual choice for solvating a single roughly-spherical solute.

In the Julia API, [`dodecahedral_unitcell`](@ref)`(T, d)` builds the same
3×3 matrix for use as the `unitcell` keyword of [`PackmolSystem`](@ref)
(paired with `unitcell_center = [cx,cy,cz]`). Since the triclinic
(parallelepiped) and dodecahedral (Wigner–Seitz) shapes are two different
choices of fundamental domain for the very same periodic cell — same
volume, same physics, just a different boundary to wrap atoms into —
[`triclinic_to_dodecahedral`](@ref) and [`dodecahedral_to_triclinic`](@ref)
convert atomic coordinates between them (analogous to GROMACS's `trjconv
-ur compact` vs `-ur rect`): useful for visualization or analysis centered
on the solute, where the more sphere-like dodecahedral shape reads more
naturally than the sheared parallelepiped.

!!! note
    Setting `pbc`/`unitcell` also implicitly confines every non-fixed structure
    type to (a slightly inflated version of) the cell, even if that structure
    declares no explicit constraint of its own. This matters because PBC by
    itself only bounds interatomic distances and how existing constraints are
    evaluated — it does not stop a molecule's own center-of-mass from drifting
    away along a direction nothing else bounds (most visibly for a half-space
    constraint like `above plane`/`below plane`, which has no periodicity of
    its own). The implicit confinement removes the need to add a redundant
    `inside box`/`inside cube` matching the cell just to keep molecules from
    wandering off.

## Box

```
inside box  xmin ymin zmin  xmax ymax zmax
outside box xmin ymin zmin  xmax ymax zmax
```

An axis-aligned rectangular box spanning the given corners.

## Cube

```
inside cube  x y z  side
outside cube x y z  side
```

An axis-aligned cube of the given `side` length, with `(x,y,z)` as its
lower corner (matching Fortran Packmol; note this differs from `box`, whose
first three numbers are also a corner but combined with a second corner
rather than a side length).

## Sphere

```
inside sphere  x y z  radius
outside sphere x y z  radius
```

A sphere of the given `radius` centered at `(x,y,z)`.

## Plane

```
above plane x y z d
over plane  x y z d
below plane x y z d
```

A half-space cut by the plane `x*a + y*b + z*c = d`, where `(a,b,c)` is the
plane's normal vector (the first three numbers) and `d` is given as the
fourth. `above`/`over` are synonyms (both accepted by the original Fortran
Packmol); `below` is the complementary half-space.

## Cylinder

```
inside cylinder  cx cy cz  vx vy vz  radius length
outside cylinder cx cy cz  vx vy vz  radius length
```

A finite, capped cylinder: `(cx,cy,cz)` is the center of one end cap,
`(vx,vy,vz)` is the axis direction (any nonzero vector — it does not need to
be normalized), `radius` is the cylinder's radius, and `length` is how far
the cylinder extends from the given end, along the axis, to the other cap.

## Ellipsoid

```
inside ellipsoid  cx cy cz  a b c  scale
outside ellipsoid cx cy cz  a b c  scale
```

An axis-aligned ellipsoid centered at `(cx,cy,cz)` with semi-axes `a`, `b`,
`c` along x, y, z respectively, uniformly scaled by the dimensionless
`scale` factor — i.e. the effective semi-axes are `a*scale`, `b*scale`,
`c*scale`. `scale` lets the same base ellipsoid (e.g. one fit to a
reference structure) be grown or shrunk without recomputing `a`, `b`, `c`.

## Rotation constraint

```
constrain_rotation <x|y|z> <center_deg> <halfwidth_deg>
```

Unlike the shape constraints above, this does not restrict where a
structure's atoms are placed — it restricts how the whole molecule is
allowed to rotate. It bounds rotation about the given axis to
`<center_deg> ± <halfwidth_deg>` degrees, as a hard bound on the optimizer's
own rotation-angle variable rather than a soft penalty. It is set once per
structure (not per constraint line, and not wrappable in `atoms ... end`);
repeat the line once per axis to constrain more than one. In the [Julia
API](julia_api.md), it is the `constrain_rotation` keyword of
[`structure_type`](@ref), taking a `Dict` such as
`Dict(:z => (0.0, 15.0))`.

## Julia constructors

Each shape above also has `Inside`/`Outside` (or `Above`/`Below`) Julia
constructors with the same parameters, for building a `PackmolSystem`
[directly from Julia code](julia_api.md) instead of a text file: `InsideBox`,
`OutsideBox`, `InsideCube`, `OutsideCube`, `InsideSphere`, `OutsideSphere`,
`InsideCylinder`, `OutsideCylinder`, `InsideEllipsoid`, `OutsideEllipsoid`,
`AbovePlane`, `BelowPlane`. Each also accepts an optional `weight` keyword
(the constraint penalty's weight in the objective function).
