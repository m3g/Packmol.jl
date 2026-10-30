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
