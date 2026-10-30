# Orthorhombic box

```
pbc a b c
pbc xmin ymin zmin  xmax ymax zmax
```

The `a b c` form is the box of the given side lengths with a corner at the
origin, from `(0, 0, 0)` to `(a, b, c)`, as in Fortran Packmol (it is the same
as `pbc 0 0 0 a b c`). The six-value form gives explicit min/max corners
instead. In both cases the center of the cell is the midpoint of the box,
`((xmin+xmax)/2, ...)`. Either way the resulting cell is diagonal (orthorhombic): no
shear between axes.
