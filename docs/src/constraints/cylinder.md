# Cylinder

```
inside cylinder  cx cy cz  vx vy vz  radius length
outside cylinder cx cy cz  vx vy vz  radius length
```

A finite, capped cylinder: `(cx,cy,cz)` is the center of one end cap,
`(vx,vy,vz)` is the axis direction (any nonzero vector — it does not need to
be normalized), `radius` is the cylinder's radius, and `length` is how far
the cylinder extends from the given end, along the axis, to the other cap.
