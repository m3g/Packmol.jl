# Constraints on some atoms

The `atom_constraints` keyword of [`structure_type`](@ref) applies constraints to some
atoms of the molecule only, as the `atoms ... end atoms` blocks of an input file. It is a
list of `atom indices => constraints` pairs (one constraint or a vector of them), with
the indices of the atoms in the structure file:

```julia
lipid = structure_type("palmitoil.pdb";
    number = 50,
    constraints = [InsideBox([0.,0.,0.],[40.,40.,14.])], # all atoms
    atom_constraints = [
        [31, 32] => BelowPlane([0.,0.,1.], 2.0), # polar head down
        [1, 2] => AbovePlane([0.,0.,1.], 12.0),  # tail up
    ],
)
```

is the input file block:

```
structure palmitoil.pdb
  number 50
  inside box 0. 0. 0. 40. 40. 14.
  atoms 31 32
    below plane 0. 0. 1. 2.
  end atoms
  atoms 1 2
    above plane 0. 0. 1. 12.
  end atoms
end structure
```

See the [Lipid bilayer](../tutorials/bilayer.md) tutorial.
