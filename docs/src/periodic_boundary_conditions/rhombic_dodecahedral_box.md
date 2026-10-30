# Rhombic dodecahedral box

`pbc dodecahedral cx cy cz d` sets up a rhombic dodecahedron cell of "size"
`d` centered at `(cx,cy,cz)` — internally just the triclinic cell `a = b =
c = d`, `α = β = 60°`, `γ = 90°` (the standard MD convention, matching
GROMACS's `editconf -bt dodecahedron -d d`), so it's written out (e.g. in
`CRYST1`) and behaves exactly like any other `unitcell`. For a given cell
volume, this shape gets a bigger minimum-image distance than a cube would
(`d` is that guaranteed minimum-image distance, in every direction) — about
70.7% of a cube's volume is enough to enclose the same sphere, which is why
it's the usual choice for solvating a single roughly-spherical solute.

In the Julia API, [`dodecahedral_unitcell`](@ref)`(T, d)` builds the same
3×3 matrix for use as the `unitcell` keyword of [`PackmolSystem`](@ref)
(paired with `unitcell_center = [cx,cy,cz]`). Since the triclinic
(parallelepiped) and dodecahedral (Wigner–Seitz) shapes are two different
choices of fundamental domain for the very same periodic cell — same
volume, same physics, just a different boundary to wrap atoms into —
[`triclinic_to_dodecahedral`](@ref) and [`dodecahedral_to_triclinic`](@ref)
convert atomic coordinates between them (analogous to GROMACS's `trjconv
-ur compact` vs `-ur rect`): useful for visualization or analysis centered
on the solute, where the more sphere-like dodecahedral shape reads more
naturally than the sheared parallelepiped.

## Example

640 water molecules (about 1 g/mL) in a rhombic dodecahedron with `d = 30` Å, shown in the dodecahedral shape:

```
pbc dodecahedral 0. 0. 0. 30.

structure water.pdb
    number 640
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl):

```@setup dodecahedral
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example dodecahedral
using Packmol, PDBTools
water = structure_type("water.pdb"; number = 640)
system = PackmolSystem([water]; tolerance = 2.0, output = "dodecahedral.pdb",
    unitcell = dodecahedral_unitcell(Float64, 30.))
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
triclinic_to_dodecahedral(system) # the shape of the output
visualize(get_atoms(system))
```
