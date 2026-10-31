# Multiple structure types

Pass one `structure_type` per kind of molecule; each gets its own
constraints:

```julia
solute = structure_type("solute.pdb"; number=1, fixed=([0.,0.,0.],[0.,0.,0.]))
water = structure_type(
    "water.pdb";
    number = 1000,
    constraints = [OutsideSphere([0.,0.,0.],15.), InsideBox([-30.,-30.,-30.],[30.,30.,30.])],
)

sys = PackmolSystem([solute, water]; output = "solvated.pdb", tolerance = 2.0)
packmol(sys)
```

## Molecules in several conformations

`structure_type` takes a vector of structure files, or a file with several models,
for several conformations of the same molecule (see [Structure
blocks](@ref "Several conformations of a molecule")): it returns then a vector of
structure types, one per conformation, among which `number` is divided as evenly as
possible. It goes into the system as a single structure type:

```julia
lipids = structure_type(["popc_1.pdb", "popc_2.pdb", "popc_3.pdb"];
    number = 60, # 20 molecules of each conformation
    atom_constraints = [[1] => AbovePlane([0., 0., 1.], 19.5), [88, 131] => BelowPlane([0., 0., 1.], 6.)],
)
water = structure_type("water.pdb"; number = 1000, constraints = [InsideBox([-30.,-30.,20.], [30.,30.,40.])])
sys = PackmolSystem([lipids, water]; output = "membrane.pdb")
```
