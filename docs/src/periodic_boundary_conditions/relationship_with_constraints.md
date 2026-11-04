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

**The centers of mass in the cell.** Minimum-image distances and explicit constraints,
by themselves, only bound how atoms interact *given* a position: neither one stops the
center of mass of a molecule from drifting along a direction that nothing else bounds. This
is most visible for a half-space constraint like `above plane`/`below plane`, which has no
periodicity of its own: a molecule above a plane could wander arbitrarily far along the
plane. So, under periodic boundary conditions, the optimizer bounds the centers of mass of
the molecules to the cell (to the box that contains it, for dodecahedral and octahedral
cells), in the packing and in the fits of the initial placement, and the centers of mass are
kept at their images in the cell. There is no need to add to every structure an `inside
box` matching the cell. The atoms of a molecule are not confined: a molecule near a face of
the cell crosses it, and interacts with the molecules on the other side through the periodic
boundary.
