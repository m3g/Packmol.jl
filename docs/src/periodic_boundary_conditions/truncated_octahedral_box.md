# Truncated octahedral box

`pbc octahedral cx cy cz d` sets up a truncated octahedron cell of "size"
`d` centered at `(cx,cy,cz)` — again internally just a triclinic cell, `a =
b = c = d`, `α = γ = acosd(1/3) ≈ 70.53°`, `β = acosd(-1/3) ≈ 109.47°`
(matching GROMACS's `editconf -bt octahedron -d d`). It's the Wigner-Seitz
cell of a BCC lattice, as opposed to the dodecahedron's FCC one, and for a
given `d` needs a bit more volume than the dodecahedron does (≈77.0% of a
cube's volume, vs. the dodecahedron's ≈70.7%) — but its shape is still
often preferred for solvating a single roughly-spherical solute.

Mirroring the dodecahedron above, [`octahedral_unitcell`](@ref)`(T, d)`
builds the matching 3×3 matrix, and [`triclinic_to_octahedral`](@ref)/
[`octahedral_to_triclinic`](@ref) convert atomic coordinates between the
skewed-parallelepiped and truncated-octahedron (Wigner-Seitz) fundamental
domains of the very same cell.
