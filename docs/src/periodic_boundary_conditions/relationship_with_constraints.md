# Relationship with constraints

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
opposite face, as if the cell tiled space. A molecule is always wrapped as a
rigid unit, never atom by atom, so it can never be torn across a periodic
boundary in a way that would leave part of it overlapping unrelated atoms.

**The shape of the cell.** For orthorhombic and general triclinic cells, the
molecules are wrapped into the parallelepiped of the cell vectors (the one
`unitcell`/CRYST1 describes), centered at `unitcell_center`. Dodecahedral and
octahedral cells (`pbc dodecahedral`/`pbc octahedral`, or a `unitcell` built with
[`dodecahedral_unitcell`](@ref)/[`octahedral_unitcell`](@ref)) are, instead, the
rhombic dodecahedron and the truncated octahedron centered at `unitcell_center`:
the molecules are wrapped to their image closest to the center, which is in that
shape. This is what the constraints of these cells refer to, and it is independent
of `periodic_boundary_style`, which only affects how `get_atoms`/`write_output`
present the *final* coordinates (see the output shapes of the [rhombic
dodecahedral](rhombic_dodecahedral_box.md#Output-shape) and [truncated
octahedral](truncated_octahedral_box.md#Output-shape) boxes). It matters for
constraints on large regions: the sphere inscribed in a truncated octahedron of size
`d` has radius `d/2`, while the one inscribed in the parallelepiped of its cell
vectors has radius `0.408 d`, so, for example, `outside sphere` constraints around a
solute (as the water around a vesicle) can use
nearly all the cell.

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
  - For a general triclinic cell it is, instead, one inward-facing
    `AbovePlane`/`BelowPlane` pair per pair of opposing cell faces (three pairs
    in total) — a box's axis-aligned edges can't represent a sheared cell, but
    three pairs of half-space cuts can.
  - For dodecahedral and octahedral cells, it is one `BelowPlane` per face of the
    rhombic dodecahedron (12) or of the truncated octahedron (14).
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
