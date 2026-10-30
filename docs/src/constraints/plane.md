# Plane

```
above plane x y z d
over plane  x y z d
below plane x y z d
```

A half-space cut by the plane `x*a + y*b + z*c = d`, where `(a,b,c)` is the
plane's normal vector (the first three numbers) and `d` is given as the
fourth. `above`/`over` are synonyms (both accepted by the original Fortran
Packmol); `below` is the complementary half-space.
