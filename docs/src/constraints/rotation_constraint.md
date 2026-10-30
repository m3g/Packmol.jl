# Rotation constraint

```
constrain_rotation <x|y|z> <center_deg> <halfwidth_deg>
```

or, with the vector of the molecule given by two sets of atoms:

```
constrain_rotation <x|y|z> <center_deg> <halfwidth_deg>
  atoms <i> <j> ...
  to <k> <l> ...
end constrain_rotation
```

Unlike the shape constraints above, this does not restrict where a
structure's atoms are placed, but how the molecule is oriented: the angle
between the given axis and a vector of the molecule must be within
`<center_deg> ± <halfwidth_deg>` degrees. The vector goes from the center of
the atoms given in the `atoms` line to the center of those given in the `to`
line (indices of the atoms in the structure file). For example, to keep the
vector from the head (atoms 1 and 2) to the tail (atom 42) of a lipid within
15° of the z axis:

```
constrain_rotation z 0. 15.
  atoms 1 2
  to 42
end constrain_rotation
```

Without the sets of atoms, the vector is the longest axis of the molecule
(its principal axis of inertia with the smallest moment, with all atoms
with the same mass). Its sign is arbitrary, so use the sets of atoms when
the direction matters (`center 90` constrains it to a plane, which doesn't
depend on the sign).

The constraint is a penalty, like the others: the squared displacement of the
tip of the vector, outside the allowed cone. Repeat the line (or block) to
constrain the angles with more than one axis. In the [Julia
API](../julia_api/index.md), it is the `constrain_rotation` keyword of
[`structure_type`](@ref).

!!! note "Difference from Fortran Packmol"
    Fortran Packmol bounds three Euler angles of the molecule, relative to
    its orientation in the structure file. The same `constrain_rotation`
    lines give different orientations here. For instance, Fortran's
    `constrain_rotation x 0 10` alone bounds the tilt of the z axis of the
    molecule (as given in its file), while here it keeps the longest axis of
    the molecule within 10° of x.
