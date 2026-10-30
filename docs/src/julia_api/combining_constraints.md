# Combining constraints

Any number of constraints can be given for a `structure_type`; a molecule
must satisfy all of them simultaneously. `InsideBox`/`OutsideSphere`/... are
regular Julia values, so they can be built programmatically (in a loop, from
a function, ...) instead of being written out by hand:

```julia
shells = [
    (r_in, r_out) -> structure_type(
        "water.pdb";
        number = 500,
        constraints = [InsideSphere([0.,0.,0.],r_out), OutsideSphere([0.,0.,0.],r_in)],
    )
    for (r_in, r_out) in [(0.,120.), (150.,200.)]
]
```
