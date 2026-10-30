# Structure blocks

```
structure <pdb file>
    number <n>
    [fixed x y z beta gamma theta]
    [connect no]
    [center | centerofmass]
    [radius <value>]
    [resnumbers <0|1|2|3>]
    [chain <c> | changechains]
    [segid <id>]
    [constrain_rotation <x|y|z> <center_deg> <halfwidth_deg>
      [atoms <i> ...
       to <j> ...
     end constrain_rotation]]
    [restart_from <file>]
    [restart_to <file>]
    <constraint lines...>
    [atoms <indices...>
        <per-atom overrides...>
    end]
end structure
```

| Keyword | Meaning |
|:--|:--|
| `number <n>` | Number of copies of this molecule to pack (required). |
| `fixed x y z beta gamma theta` | Don't pack this molecule — place it once at this position and (Euler angle) orientation and leave it fixed. `number` is optional for a fixed structure (it is a single molecule). See [Fixed molecules](../constraints/fixed.md). |
| `connect no` | Don't write the connectivity (`CONECT` records) of this structure to the output. |
| `center` / `centerofmass` | Only valid together with `fixed`: re-center the molecule (geometric center or mass-weighted center of mass) before applying the fixed position/orientation. Without them, the molecule is rotated about the origin of its coordinates, and then translated. (Fortran Packmol uses unit masses for `centerofmass`, so its center of mass is the geometric center.) |
| `radius <value>` | Atom radius for every atom of this structure, overriding the default of `tolerance/2`. |
| `resnumbers <n>` | Residue numbering scheme for the output — see [Residue numbering](@ref) below. |
| `chain <c>` | Fix every molecule of this type to chain identifier `<c>` in the output, instead of the default auto-assigned one. |
| `changechains` | Alternate between two auto-assigned chain identifiers every other molecule, instead of one chain for the whole type. Not compatible with `chain`. |
| `segid <id>` | Segment identifier (at most 4 characters) for every atom of this type, in columns 73-76 of the output PDB file. Without `segid`, the segment identifiers of the structure's own PDB file are kept (blank if it has none). |
| `constrain_rotation <axis> <center_deg> <halfwidth_deg>` | Keep the angle between `<axis>` (`x`, `y` or `z`) and a vector of the molecule within `<center_deg> ± <halfwidth_deg>` degrees: the longest axis of the molecule, or, followed by `atoms <i> ...`, `to <j> ...` and `end constrain_rotation` lines, the vector between the centers of two sets of atoms. See [Rotation constraint](@ref). |
| `restart_from <file>` | Restart only this structure type's molecules — see [Restarting a run](@ref). |
| `restart_to <file>` | Write the coordinates of this structure type's molecules only — see [Restarting a run](@ref). |
| `atoms <indices> ... end` | Restrict the enclosed lines (`radius`, constraints) to the listed 1-based atom indices instead of the whole molecule. |

Constraint lines (`inside box`, `outside sphere`, `above plane`, ...) are
documented on the [Constraints](../constraints/index.md) page — any number of them can
appear in a structure block, optionally inside an `atoms ... end` block.
