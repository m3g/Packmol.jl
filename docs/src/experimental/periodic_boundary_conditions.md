# Periodic boundary conditions

!!! warning "Experimental"
    This documents the native Julia packing engine's periodic boundary
    syntax, shared between [input files](input_files.md) and the
    [Julia API](julia_api.md). It tracks the original Fortran
    [Packmol](http://github.com/m3g/packmol) keywords where they exist, but
    is not yet a complete port — see [Home](index.md) for what's covered so
    far.

```
pbc a b c
pbc xmin ymin zmin  xmax ymax zmax
pbc dodecahedral cx cy cz d
pbc octahedral cx cy cz d
unitcell a b c alpha beta gamma
```

`pbc`/`unitcell` are not shape constraints declared inside a `structure`
block (those are documented on the [Constraints](constraints.md) page) —
they are global keywords, given once per input file (or once per
`PackmolSystem`, via the `unitcell`/`unitcell_center` keywords in the [Julia
API](julia_api.md)), that set up a periodic simulation cell. Every form
above ends up stored the same way: as a 3×3 `unitcell` matrix (its columns
are the cell vectors) plus a `unitcell_center`. `pbc` is a convenience for
the two most common shapes (an orthorhombic box, or a rhombic dodecahedron/
truncated octahedron computed from a single size parameter); `unitcell`
gives the fully general triclinic cell directly.

## Orthorhombic box

```
pbc a b c
pbc xmin ymin zmin  xmax ymax zmax
```

The `a b c` form is a box of the given side lengths, centered at the
origin. The six-value form gives explicit min/max corners instead — the box
is then centered at their midpoint, `((xmin+xmax)/2, ...)`, rather than at
the origin. Either way the resulting cell is diagonal (orthorhombic): no
shear between axes.

## Triclinic box

```
unitcell a b c alpha beta gamma
```

The fully general periodic cell, CRYST1-style: side lengths `a`, `b`, `c`
and angles `alpha`, `beta`, `gamma` (in degrees, between the corresponding
pairs of sides), built with the PDB convention (`a` along x, `b` in the
xy-plane). Always centered at the origin from the input file; the Julia API
accepts any `unitcell_center` alongside an explicit 3×3 `unitcell` matrix
built however the caller likes (`dodecahedral_unitcell`/`octahedral_unitcell`
below build exactly this kind of matrix from a single size parameter, as a
convenience).

An orthorhombic `pbc` is just the special case `alpha = beta = gamma = 90°`
— everything below (minimum-image distances, implicit confinement, CRYST1
output) is expressed in terms of the general triclinic cell and reduces to
the simpler orthorhombic behavior automatically when the cell happens to be
diagonal.

## Rhombic dodecahedral box

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

## Truncated octahedral box

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

## Choosing the output shape

Called on a `PackmolSystem` directly (no coordinates involved),
`triclinic_to_dodecahedral(packmol_system)`/`triclinic_to_octahedral(packmol_system)`
(and their `..._to_triclinic` inverses) set
`packmol_system.periodic_boundary_style` (`:dodecahedral`/`:octahedral`, or
`:triclinic`, the default) and return `packmol_system` — the shape that
[`get_atoms`](@ref) and [`write_output`](@ref) wrap molecule positions into
when they materialize atomic coordinates, generated fresh from
`packmol_system.molecule_positions` each time rather than stored. No
coordinates are touched by the toggle itself; only the next
`get_atoms`/`write_output` call is affected, and only that call's *output*
— see [Relationship with constraints](@ref) below for why packing itself is
unaffected by this choice.

## Relationship with constraints

Once `pbc`/`unitcell` is set, it changes how every distance in the system
is computed, and how every explicit constraint (declared inside a
`structure ... end structure` block, or via `Inside`/`Outside` constructors
in the Julia API) is evaluated — on top of adding one implicit constraint of
its own.

**Interatomic distances.** With PBC active, the tolerance/overlap distance
between two atoms is always the *minimum-image* distance: each molecule's
center of mass is wrapped to the periodic image nearest the cell's
center (`unitcell_center`) before distances are computed, so a molecule
near one face of the cell interacts naturally with molecules near the
opposite face, as if the cell tiled space. This wrapping always uses the
ordinary skewed-parallelepiped image (the same one `unitcell`/CRYST1
describes) — it is entirely independent of `periodic_boundary_style`, which
only affects how `get_atoms`/`write_output` present the *final* coordinates
(see above). A molecule is always wrapped as a rigid unit, never atom by
atom, so it can never be torn across a periodic boundary in a way that
would leave part of it overlapping unrelated atoms.

**Explicit constraints.** Any `inside`/`outside`/`above`/`below` constraint
declared for a structure is checked against that same wrapped position —
so, for example, `inside sphere cx cy cz r` together with PBC confines a
molecule to a sphere drawn in the *wrapped* cell, not in unwrapped/raw
coordinate space.

**Implicit confinement (the boundary added).** Setting `pbc`/`unitcell`
also implicitly confines every non-fixed structure type to (a slightly
inflated version of) the cell, even if that structure declares no explicit
constraint of its own:

  - For an orthorhombic cell this is a single `InsideBox` spanning the
    cell.
  - For a general triclinic cell (including the dodecahedral and octahedral
    ones above) it is, instead, one inward-facing `AbovePlane`/`BelowPlane`
    pair per pair of opposing cell faces (three pairs in total) — a box's
    axis-aligned edges can't represent a sheared cell, but three pairs of
    half-space cuts can.
  - Each face is pushed outward by `tolerance / 2`, so the confining
    boundary ends up `tolerance` larger than the literal cell, in total,
    along each axis: an atom sitting exactly on the true cell boundary
    still needs its own radius' worth of room (the default atom radius is
    `tolerance / 2`) to interact naturally with its periodic image across
    that boundary, the same way it would with a same-radius neighbor
    anywhere in the interior. Without this slack, the implicit confinement
    would push back *before* that natural interaction has a chance to
    happen, visibly flattening the packing right at the cell boundary
    instead of letting it look like an unbroken, continuous interior.

This implicit confinement exists because minimum-image distances and
explicit constraints, by themselves, only bound how atoms interact *given*
a position — neither one stops a molecule's own center of mass from
drifting away indefinitely along a direction nothing else bounds. This is
most visible for a half-space constraint like `above plane`/`below plane`,
which has no periodicity of its own: without the implicit confinement,
nothing would stop a molecule from wandering arbitrarily far along the
plane's unconstrained directions, or off to its own unconstrained side. The
implicit confinement removes the need to add a redundant `inside
box`/`inside cube` matching the cell to every structure just to keep
molecules from wandering off; it applies automatically to every free
(non-`fixed`) structure type, on top of whatever explicit constraints that
structure already declares.
