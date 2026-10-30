# Julia constructors

Each shape above also has `Inside`/`Outside` (or `Above`/`Below`) Julia
constructors with the same parameters, for building a `PackmolSystem`
[directly from Julia code](../julia_api/index.md) instead of a text file: `InsideBox`,
`OutsideBox`, `InsideCube`, `OutsideCube`, `InsideSphere`, `OutsideSphere`,
`InsideCylinder`, `OutsideCylinder`, `InsideEllipsoid`, `OutsideEllipsoid`,
`AbovePlane`, `BelowPlane`, `AboveGaussian`, `BelowGaussian`, `AboveSin`,
`BelowSin`, `AboveCos`, `BelowCos`, `AboveExponential`, `BelowExponential`,
`AboveRadialGaussian`, `BelowRadialGaussian`, `AboveRadialSin`,
`BelowRadialSin`, `AboveRadialCos`, `BelowRadialCos`,
`AboveRadialExponential`, `BelowRadialExponential`. Each also accepts an
optional `weight` keyword (the constraint penalty's weight in the objective
function).

For boxes and cubes, the positional arguments are those of the input file:
`InsideBox([xmin, ymin, zmin], [xmax, ymax, zmax])` and
`InsideCube([xmin, ymin, zmin], side)`. They can also be given by their center, with
keywords: `InsideBox(center=[x, y, z], sides=[a, b, c])` and
`InsideCube(center=[x, y, z], side=a)`.
