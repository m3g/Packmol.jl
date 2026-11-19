# Packing by compression

Some systems do not pack from a random placement at their density: long chains,
polymers, crowded large molecules. The molecules jam before they fill the space, each
blocked by its neighbors, while the same molecules pack in seconds at a lower density.
The keyword `compress` (Packmol.jl only) packs such systems at a lower density and then
compresses them, in steps, to the density of the input:

```
compress              # pack at half the density, then compress
```

The input describes the system at its density (the constraints, the periodic cell, the
numbers of molecules), and `compress [<fraction> [<step>]]` takes two optional values:

- `fraction` (default `0.5`): the fraction of that density at which the system is first
  packed. With `0.5`, the system is packed in twice its volume, at half the density; with
  `0.3`, at 30% of it.
- `step` (default `0.05`): the initial step of the compression, the fraction by which
  the linear dimensions shrink at each step (5%).

`compress` alone is `compress 0.5 0.05`:

1. The constraints and the periodic cell are expanded by `fraction^(-1/3)` (by 1.26 for
   `0.5`), about the center of the cell (the origin, without periodic boundary
   conditions), and the system is packed there.
2. At each step, the geometry is scaled by `1 - step` about the same center: the
   constraints, the cell, and the positions of the free molecules, each moved as a whole
   (the center of mass of a molecule is scaled; a polymer moves with the center of its
   monomers). Fixed molecules stay where they are. The system is packed again, from there, for up to 100 loops.
3. A step that does not converge is undone, and tried again with half the step; after
   two steps that converge, the step grows by 1.5, up to the initial step. The
   compression ends at the geometry of the input, or stops (with a warning, the output
   at the density reached) if a step smaller than 0.0005 does not converge.

The report shows each step:

```
  Compression step    1: density   58.32% of the target, not converged (4.3 s)
  Compression step    2: density   53.95% of the target, converged (0.2 s)
  ...
  Compressed to the density of the input in 20 steps (41.3 s).
```

The compression works with rigid molecules and with [Polymers](@ref). With polymers
it goes much further: the chains of the condensate of the tutorial
[Biomolecular condensate](@ref), packed as rigid molecules, jam at about 220 mg/mL,
while released as flexible chains they are compressed to 400 mg/mL in about 40 seconds.
Small steps matter: compressed in steps of 3%, without the adaptive step, they stopped at
263 mg/mL (the chains pushed into each other by a large step get threaded through each
other's monomers), and in steps of 0.5% they reached 400 mg/mL in 70 steps.

In Julia, `compress` is a keyword argument of [`packmol`](@ref) and of
[`PackmolSystem`](@ref): `compress = true` (the defaults), `compress = 0.3` (the
fraction), or `compress = (0.3, 0.02)` (the fraction and the initial step).

## Periodic cells and fixed molecules

All periodic cells are supported (orthorhombic, triclinic, dodecahedral, octahedral):
they are scaled uniformly about their center, which keeps their shape.

Fixed molecules are supported, and are not moved: the free molecules are packed around
them in the expanded geometry and compressed toward them. This is right when the free
molecules are placed freely around the fixed ones (polymers or a solvent around a
nanoparticle or a protein), since the expanded cell and constraints contain them.

!!! note "Limitations"
    - The constraints must be boxes, cubes, spheres, cylinders, ellipsoids or planes (the
      others cannot be scaled).
    - The constraints are scaled and the fixed molecules are not: a constraint that is
      placed relative to a fixed molecule loses its relation to it in the expanded
      geometry, and recovers it only at the end of the compression. For example, the
      membrane of the [Zika virion (a large system)](@ref), in spherical shells chosen
      to fit inside the fixed shell of proteins, crossed by their transmembrane segments:
      expanded, the shells of the membrane would be beyond the proteins, and the lipids
      packed there would be compressed through them. In such systems, pack
      the free molecules whose regions are tied to the fixed ones without `compress`.
    - Not with `restart_from`.
