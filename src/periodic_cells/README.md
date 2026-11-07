# periodic_cells

Periodic boundary conditions. The cell is a matrix whose columns are the cell vectors,
centered at `unitcell_center`.

| File | Contents |
|---|---|
| `periodic_cells.jl` | The shared functions of the periodic cells: the nearest periodic image, the check that fixed molecules are inside the cell, the Wigner-Seitz cells (rhombic dodecahedron and truncated octahedron, detected from the cell matrix), the half-extent of the cell that bounds the centers of mass, the wrapping of the stored centers of mass after each loop, and `unitcell_matrix` (the cell from a b c α β γ). |
| `dodecahedron.jl` | The rhombic dodecahedron: `dodecahedral_unitcell`, and the conversions between the triclinic cell and the shape of the dodecahedron. |
| `octahedron.jl` | The truncated octahedron: `octahedral_unitcell`, and its conversions. |
| `wrapping.jl` | Wrapping of positions (`wrap_to_center`) and of centers of mass (`_wrap_cm`, to the image in the Wigner-Seitz cell) into the cell. |

The constraints are evaluated at the image of each molecule in the cell, and the centers
of mass are bounded to the cell in the optimizations (see `../initial_approximation/` and
`../packing/`).
