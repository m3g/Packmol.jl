# Packing by compression

Some systems do not pack from a random placement at their density: long chains,
polymers, crowded large molecules. The molecules jam before they fill the space, each
blocked by its neighbors, while the same molecules pack in seconds at a lower density.
The keyword `use_compression` packs such systems at a lower density, and compresses
them from there, in steps, to the density of the input:

```
use_compression       # pack at half the density, then compress
```

The input describes the system at its density (the constraints, the periodic cell, the
numbers of molecules), and `use_compression [<fraction> [<step>]]` takes two optional
values:

- `fraction` (default `0.5`): the fraction of that density at which the system is first
  packed. With `0.5`, the system is packed in twice its volume, at half the density; with
  `0.3`, at 30% of it.
- `step` (default `0.05`): the initial step of the compression, the fraction by which
  the linear dimensions shrink at each step (5%).

`use_compression` alone is `use_compression 0.5 0.05`:

1. The constraints and the periodic cell are expanded by `fraction^(-1/3)` (by 1.26 for
   `0.5`), about the center of the cell (the origin, without periodic boundary
   conditions), and the system is packed there.
   If the packing at that density does not converge, the geometry is expanded further
   (to half the density, the molecules with it) and the system packed again from there,
   up to three times.
2. At each step, the geometry is scaled by `1 - step` about the same center: the
   constraints, the cell, and the positions of the free molecules, each moved as a whole
   (the center of mass of a molecule is scaled; a flexible molecule, a polymer or a
   molecule with rotatable bonds, moves with the center of its fragments), and the fixed
   molecules, scaled with the radii of their atoms (see below). The system is packed again, from there, for up to 100 loops.
3. When a step does not converge, the system breathes: the packing of the failed step is
   expanded to one step below the last density that converged, and packed there (a
   relaxation), keeping the rearrangements the failed step made. The compression
   resumes from the relaxed packing, with half the step (if the relaxation does not
   converge either, or after 20 relaxations, it resumes from the last converged packing).
   After two steps that converge, the step grows by 1.5, up to the initial step. The
   compression ends at the geometry of the input, or stops (with a warning, the output
   at the density reached) if a step smaller than 0.0005 does not converge.

The report shows each step:

```
  Compression step    1: density   58.32% of the target, converged (0.8 s)
  Compression step    2: density   68.02% of the target, converged (1.5 s)
  Compression step    3: density   79.33% of the target, converged (3.0 s)
  Compression step    4: density   92.53% of the target, not converged (7.6 s)
  Relaxation    1:        density   68.02% of the target, converged (0.1 s)
  Compression step    5: density   73.39% of the target, converged (0.1 s)
  ...
  Compression step   15: density  100.00% of the target, converged (1.9 s)

  Compressed to the density of the input in 15 steps (28.0 s).
```

Undoing a failed step and halving it, without relaxing it, stopped the compression of the
particle of the [Lignin particle in acetone and water](@ref) tutorial short of the
density in most runs: the next steps failed at the same configuration. Relaxing at
every failed step, it reached the density in all runs (5 of 5 seeds), while relaxing
only as a last resort, once the step had shrunk below its minimum, failed in one of five
and took 1.6 times longer.

The compression works with rigid molecules and with flexible ones: [Polymers](@ref) and
molecules with [Rotatable bonds](@ref). With flexible molecules it goes much further, as
they change their shape to fill the space: the chains of the condensate of the tutorial
[Biomolecular condensate](@ref), packed as rigid molecules, jam at about 220 mg/mL,
while as flexible chains they are compressed to 400 mg/mL in about a minute; the lignin
oligomers of the [Lignin particle in acetone and water](@ref), flexible at every
rotatable bond (`rotatable all`), reach the density of lignin, 1.30 g/cm³. Large steps
push the flexible molecules into each other, and leave them entangled (chains threaded
through the fragments of others), hence the adaptive step and the relaxations.

In Julia, `use_compression` is a keyword argument of [`packmol`](@ref) and of
[`PackmolSystem`](@ref): `use_compression = true` (the defaults), `use_compression = 0.3` (the
fraction), or `use_compression = (0.3, 0.02)` (the fraction and the initial step).

## Periodic cells and fixed molecules

All periodic cells are supported (orthorhombic, triclinic, dodecahedral, octahedral):
they are scaled uniformly about their center, which keeps their shape.

Fixed molecules are supported, and are scaled with the geometry: the coordinates of their
atoms about the center of the scaling, and the radii of their atoms, such that the
constraints placed relative to them keep their relation to them. At the end of the
compression, they are back at the coordinates and with the radii of the input. For
example, the [Zika virion (a large system)](@ref), packed with `use_compression 0.5`: its
membrane, in spherical shells chosen to fit inside the fixed shell of proteins, crossed by
their transmembrane segments, swells and shrinks with the proteins, and the virion is
compressed to its density in 5 steps (in about 2 minutes, while it packs without
compression in about one). With the proteins kept at their place, the shells of the
membrane, expanded, were beyond them, and the compression had not gone beyond a third of
the density after 25 minutes.

!!! note "Fixed molecules, scaled"
    A fixed molecule is scaled as a whole: its atoms move apart, and grow, by the same
    factor (by 1.26 for `use_compression 0.5`). The free molecules around it are not, so in
    the expanded geometry they meet a larger version of it: free molecules placed in its
    cavities or grooves (water inside a protein, ligands in a pocket, the lipids crossed by
    a transmembrane helix) see them as they are only at the end of the compression, and may
    need the last steps to settle. This works when the free molecules are around the fixed
    ones, or in regions whose shapes follow them; it is not known to work in every case. If
    the compression stops short of the density, the output is written at the density
    reached, the fixed molecules scaled with it.

!!! note "Limitations"
    - The constraints must be boxes, cubes, spheres, cylinders, ellipsoids or planes (the
      others cannot be scaled).
    - Not with `restart_from`.
