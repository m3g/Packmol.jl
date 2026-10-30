# Gaussian, wave, and exponential surfaces

```
above gaussian    ux uy uz  ax ay az  d0 amplitude center sigma
below gaussian    ux uy uz  ax ay az  d0 amplitude center sigma
above sin         ux uy uz  ax ay az  d0 amplitude wavelength phase
below sin         ux uy uz  ax ay az  d0 amplitude wavelength phase
above cos         ux uy uz  ax ay az  d0 amplitude wavelength phase
below cos         ux uy uz  ax ay az  d0 amplitude wavelength phase
above exponential ux uy uz  ax ay az  d0 amplitude center rate
below exponential ux uy uz  ax ay az  d0 amplitude center rate
```

These generalize `plane` into a curved half-space `up·x = h(along·x)`: a
Gaussian bump, sine/cosine wave, or exponential ramp, extruded along the
direction perpendicular to `along` (within the plane orthogonal to `up`).
`(ux,uy,uz)` (`up`) plays the same role as `plane`'s normal and is used
un-normalized, exactly as there; `(ax,ay,az)` (`along`) picks the direction
the profile varies along and is normalized internally, so it need not be
given as a unit vector. `d0` is the base offset (as `plane`'s `d`); with
`amplitude = 0` each of these degenerates exactly to `plane`. `cos` is `sin`
with its `phase` shifted by `pi/2` internally — there is no separate
underlying shape.

- **gaussian**: `h(s) = d0 + amplitude * exp(-(s - center)^2 / (2*sigma^2))`
- **sin** / **cos**: `h(s) = d0 + amplitude * sin(2*pi/wavelength * s + phase)`
  (`cos` adds `pi/2` to `phase`)
- **exponential**: `h(s) = d0 + amplitude * exp(rate * (s - center))`

where `s = along·x` (using the normalized `along`).


## Example

390 water molecules (about 1 g/mL) in a box, above a Gaussian bump of height 6 Å and width 4 Å along x:

```
structure water.pdb
    number 390
    inside box -15. -15. 0. 15. 15. 15.
    above gaussian 0. 0. 1. 1. 0. 0. 0. 6. 0. 4.
end structure
```

The same system, packed with the [Julia API](../julia_api/index.md) and shown with the
`visualize` function of [PDBTools.jl](https://github.com/m3g/PDBTools.jl):

```@setup gaussian
using Packmol
cp(joinpath(Packmol.RecipesDirectory, "test", "data", "water.pdb"), "water.pdb"; force = true)
```

```@example gaussian
using Packmol, PDBTools
water = structure_type("water.pdb";
    number = 390,
    constraints = [InsideBox([-15., -15., 0.], [15., 15., 15.]), AboveGaussian([0., 0., 1.], [1., 0., 0.], 0., 6., 0., 4.)],
)
system = PackmolSystem([water]; tolerance = 2.0, output = "gaussian.pdb")
redirect_stdio(stdout = devnull, stderr = devnull) do # hide
packmol(system)
end # hide
visualize(get_atoms(system); unitcell = UnitCell([30., 30., 15.]), unitcell_origin = [-15., -15., 0.])
```

## Radial variants

```
above radial_gaussian    ux uy uz  cx cy cz  d0 amplitude sigma
below radial_gaussian    ux uy uz  cx cy cz  d0 amplitude sigma
above radial_sin         ux uy uz  cx cy cz  d0 amplitude wavelength phase
below radial_sin         ux uy uz  cx cy cz  d0 amplitude wavelength phase
above radial_cos         ux uy uz  cx cy cz  d0 amplitude wavelength phase
below radial_cos         ux uy uz  cx cy cz  d0 amplitude wavelength phase
above radial_exponential ux uy uz  cx cy cz  d0 amplitude rate
below radial_exponential ux uy uz  cx cy cz  d0 amplitude rate
```

Radially symmetric counterparts of the three shapes above: the surface is
`up·x = h(r)`, where `r` is the radial distance from the axis line through
`center` (`(cx,cy,cz)`) parallel to `up`, measured perpendicular to `up`
(exactly as `cylinder`'s radial distance from its axis). `up` is normalized
internally here, since it doubles as that axis. `radial_gaussian` gives a
smooth dome; `radial_sin`/`radial_cos` give concentric ripples; and
`radial_exponential` gives a cone-like spike — these last two have a genuine
kink at `r = 0` (their gradient's radial component is taken to be zero
exactly on the axis).

- **radial_gaussian**: `h(r) = d0 + amplitude * exp(-r^2 / (2*sigma^2))`
- **radial_sin** / **radial_cos**: `h(r) = d0 + amplitude * sin(2*pi/wavelength * r + phase)`
  (`radial_cos` adds `pi/2` to `phase`)
- **radial_exponential**: `h(r) = d0 + amplitude * exp(rate * r)`
