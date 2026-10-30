# Rotation constraints

[`structure_type`](@ref)'s `constrain_rotation` keyword keeps the angle between
an axis and a vector of the molecule within `center ± halfwidth` degrees (see
[Rotation constraint](@ref)). Each axis takes `(center, halfwidth)`, for the
longest axis of the molecule, or `(center, halfwidth, from_atoms, to_atoms)`,
for the vector between the centers of two sets of atoms:

```julia
lipid = structure_type(
    "lipid.pdb";
    number = 100,
    # head (atoms 1 and 2) to tail (atom 42) within 15° of z
    constrain_rotation = Dict(:z => (0.0, 15.0, [1, 2], [42])),
    constraints = [InsideBox([0.,0.,0.],[80.,80.,40.])],
)
```
