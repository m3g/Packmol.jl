# Placement on a grid (Bentopy)

!!! note "A port of Bentopy's algorithm"
    The function on this page is a direct port to Julia of the placement algorithm of
    **[Bentopy](https://github.com/marrink-lab/bentopy)** (its packing step,
    `bentopy-pack`), by Westendorp, Stevens, Brown, Dommer, Wassenaar, Bruininks and
    Marrink:

    > M. S. S. Westendorp et al., *Compartment-guided assembly of large-scale molecular
    > models with Bentopy*, Protein Science (2026) e70480.
    > [doi:10.1002/pro.70480](https://doi.org/10.1002/pro.70480)

    Bentopy builds models of whole cells, tens of thousands of proteins, RNAs and
    metabolites, in minutes, by placing the molecules one after the other on an occupancy
    grid, with no optimization at all. Packing its published model of the JCVI-Syn3A
    minimal cell is what showed us how well this works as a starting point for Packmol,
    and what led to `voxel_place`. The algorithm, and its effectiveness, are Bentopy's:
    if you use `voxel_place`, please cite Bentopy.

## What it does

[`voxel_place`](@ref), the Julia port of Bentopy's packing algorithm, takes a system (or
an input file), places its free molecules as Bentopy does, and drops the molecules that
do not fit:

- The space is a grid of voxels (of side `resolution`): the periodic cell or, without
  periodic boundary conditions, the box of the regions of the molecules. The voxels
  taken by the fixed molecules are marked.
- Each molecule is voxelized conservatively at a random rotation: the voxels that
  intersect a sphere of the radius of each of its atoms (Bentopy's bead radius). Two
  molecules that do not share voxels are therefore at distances larger than the sum of
  the radii of their atoms: the tolerance.
- The compartment of each type (Bentopy's masks) is the set of voxels whose corners
  satisfy the constraints of the type.
- The types are placed in decreasing order of their moment (the sum of the squared
  distances of the atoms to their center), each type completely before the next.
- Each try puts the molecule at a random free voxel of its compartment, and accepts it
  if none of its voxels is taken. The rotation changes after each placement, and after
  `max_tries_mult × number ÷ max_tries_rot_div` unsuccessful tries; the type gives up
  after `max_tries_mult × number` unsuccessful tries (defaults 1000 and 100, as in
  Bentopy), and the molecules that were not placed are **dropped**.

```julia
using Packmol
placed, dropped = voxel_place("cell.inp")   # or voxel_place(system)
packmol(placed)
```

In an input file, the keyword `voxel_place` (optionally followed by `max_tries_mult`) does
the same within `packmol`, also from the command line: the free molecules are placed by
`voxel_place`, and the system of those that fit is packed.

```
tolerance 4.0
voxel_place 100000
...
```

The molecules that do not fit are reported with a warning, and with a **significant
warning** when more than 10% of the molecules of a structure type are left out. The
voxel placement puts rigid molecules at random positions and rotations, and can fall far
short of the space available: for the lipids of the Zika virion membrane (heads and tails
constrained to opposite sides of a surface), it placed a quarter to a half of the lipids
that the volume of the leaflets (`free_volume`) could hold. Compare the numbers with the
volume available before relying on the placement, and pack without it in such cases.

`placed` is a system with the molecules that fit, which starts from their positions
(`placed.restart_from`, see [Restart files](@ref)), and `dropped` is the number of molecules
of each structure type that did not fit (the types are listed at the end of the placement).
The composition of `placed` is not that of the input when some molecules do not fit: check
it before using the result, and enlarge the region or reduce the numbers if the difference
matters.

The `resolution` is by default 2.5 times the largest radius of the atoms of the free
molecules (Bentopy's models use voxels of 5 Å and beads of 2 Å radius: a Packmol
tolerance of 4 Å). A `radius` keyword sets the same radius for all atoms, as Bentopy's
bead radius.

The port takes the regions of the molecules from Packmol's constraints instead of
Bentopy's compartment files, which asks for a few adaptations: constraints that apply to some
atoms of a molecule only, and orientation constraints, are checked on the atoms of a
try that does not collide; without periodic boundary conditions each type needs a
constraint that bounds its region (inside a box, cube, sphere, ellipsoid or cylinder);
periodic cells must be orthorhombic. It works in 3D only.

## When to use it

For crowded systems of many large molecules, where the usual initial placement leaves
large molecules interpenetrated and the packing stalls with a few of them that cannot
find room. On Bentopy's model of the JCVI-Syn3A minimal cell (the 41,967 proteins, RNAs
and metabolites of Bentopy's input, 16 million coarse-grained beads, inside a sphere of
950 Å, with the membrane and the chromosome fixed, tolerance 4 Å; see the
[Minimal cell](@ref "Minimal cell (voxel placement of Bentopy)") tutorial):

- the usual initial placement did not converge in 40 loops;
- `voxel_place` (with Bentopy's `max_tries_mult = 100_000` of this model) placed 41,835
  molecules in 16 seconds, leaving out 132 RNAs, as Bentopy did, and the packing from it
  converged at the first loop. Solvating it then (17.3 million water beads and 3 million
  glycerol placeholders, the cell fixed, tolerance 2 Å) took 2 loops: the whole cell
  built in about 11 minutes.

On a part of the same cell (3,122 molecules in a sphere of 400 Å), Bentopy and
`voxel_place` placed the same molecules (both dropped the same RNA) in a few seconds,
and the packings from both placements converged at the first loop.

For solvents and dilute solutions it is not needed: the usual initial placement and the
packing are fast there.

## Bentopy and Packmol together

Bentopy itself is the tool to build models of cells and other crowded compartments:
compartments from voxel masks and combinations of them, concentrations, and the
features its authors keep adding. Check [its repository](https://github.com/marrink-lab/bentopy)
for what it offers before building such a system by other means. Its output can be
taken by Packmol, for what Packmol adds:

- **Removing the remaining overlaps** at a given tolerance: restart from Bentopy's
  coordinates (`restart_from`, a coordinate file with the atoms in the order of the
  structure blocks of the input, each structure block a molecule type of Bentopy's
  model). On the minimal cell this took 2 loops (94 s).
- **Solvation**, with any solvent and mixtures, the molecules of the model fixed: the
  minimal cell was solvated with 17.3 million water beads and 3 million glycerol
  placeholders in 5.4 minutes at a tolerance of 2 Å.
- **Constraints** that Bentopy does not express: regions of any shape (planes, curved
  surfaces, surfaces given by points), constraints on some atoms of a molecule, and
  orientations.

Bentopy, like `voxel_place`, leaves out the molecules that do not fit (on the minimal
cell, 132 of 686 RNAs). Whether Packmol can add them to the placed system, the others
free to move, is an open question.
