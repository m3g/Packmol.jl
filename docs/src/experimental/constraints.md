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

!!! note
    Constraints are independent of, but interact with, periodic boundary
    conditions: setting `pbc`/`unitcell` changes how every constraint below
    is evaluated (atoms are wrapped to the periodic image nearest the cell
    center first), and adds one implicit confinement constraint of its own
    on top of whatever's declared explicitly. See [Periodic boundary
    conditions](periodic_boundary_conditions.md) for the `pbc`/`unitcell`
    keywords themselves and the full relationship between the two.

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

## Gaussian, wave, and exponential surfaces

```
above gaussian    ux uy uz  ax ay az  d0 amplitude center sigma
below gaussian    ux uy uz  ax ay az  d0 amplitude center sigma
above sin         ux uy uz  ax ay az  d0 amplitude wavelength phase
below sin         ux uy uz  ax ay az  d0 amplitude wavelength phase
above cos         ux uy uz  ax ay az  d0 amplitude wavelength phase
below cos         ux uy uz  ax ay az  d0 amplitude wavelength phase
above exponential ux uy uz  ax ay az  d0 amplitude center rate
below exponential ux uy uz  ax ay az  d0 amplitude center rate
```

These generalize `plane` into a curved half-space `up·x = h(along·x)`: a
Gaussian bump, sine/cosine wave, or exponential ramp, extruded along the
direction perpendicular to `along` (within the plane orthogonal to `up`).
`(ux,uy,uz)` (`up`) plays the same role as `plane`'s normal and is used
un-normalized, exactly as there; `(ax,ay,az)` (`along`) picks the direction
the profile varies along and is normalized internally, so it need not be
given as a unit vector. `d0` is the base offset (as `plane`'s `d`); with
`amplitude = 0` each of these degenerates exactly to `plane`. `cos` is `sin`
with its `phase` shifted by `pi/2` internally — there is no separate
underlying shape.

- **gaussian**: `h(s) = d0 + amplitude * exp(-(s - center)^2 / (2*sigma^2))`
- **sin** / **cos**: `h(s) = d0 + amplitude * sin(2*pi/wavelength * s + phase)`
  (`cos` adds `pi/2` to `phase`)
- **exponential**: `h(s) = d0 + amplitude * exp(rate * (s - center))`

where `s = along·x` (using the normalized `along`).

### Radial variants

```
above radial_gaussian    ux uy uz  cx cy cz  d0 amplitude sigma
below radial_gaussian    ux uy uz  cx cy cz  d0 amplitude sigma
above radial_sin         ux uy uz  cx cy cz  d0 amplitude wavelength phase
below radial_sin         ux uy uz  cx cy cz  d0 amplitude wavelength phase
above radial_cos         ux uy uz  cx cy cz  d0 amplitude wavelength phase
below radial_cos         ux uy uz  cx cy cz  d0 amplitude wavelength phase
above radial_exponential ux uy uz  cx cy cz  d0 amplitude rate
below radial_exponential ux uy uz  cx cy cz  d0 amplitude rate
```

Radially symmetric counterparts of the three shapes above: the surface is
`up·x = h(r)`, where `r` is the radial distance from the axis line through
`center` (`(cx,cy,cz)`) parallel to `up`, measured perpendicular to `up`
(exactly as `cylinder`'s radial distance from its axis). `up` is normalized
internally here, since it doubles as that axis. `radial_gaussian` gives a
smooth dome; `radial_sin`/`radial_cos` give concentric ripples; and
`radial_exponential` gives a cone-like spike — these last two have a genuine
kink at `r = 0` (their gradient's radial component is taken to be zero
exactly on the axis).

- **radial_gaussian**: `h(r) = d0 + amplitude * exp(-r^2 / (2*sigma^2))`
- **radial_sin** / **radial_cos**: `h(r) = d0 + amplitude * sin(2*pi/wavelength * r + phase)`
  (`radial_cos` adds `pi/2` to `phase`)
- **radial_exponential**: `h(r) = d0 + amplitude * exp(rate * r)`

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
`AbovePlane`, `BelowPlane`, `AboveGaussian`, `BelowGaussian`, `AboveSin`,
`BelowSin`, `AboveCos`, `BelowCos`, `AboveExponential`, `BelowExponential`,
`AboveRadialGaussian`, `BelowRadialGaussian`, `AboveRadialSin`,
`BelowRadialSin`, `AboveRadialCos`, `BelowRadialCos`,
`AboveRadialExponential`, `BelowRadialExponential`. Each also accepts an
optional `weight` keyword (the constraint penalty's weight in the objective
function).
