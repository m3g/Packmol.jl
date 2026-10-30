# Tutorials

The first tutorials are the classic examples of
[Packmol](https://m3g.github.io/packmol/examples.shtml), followed by systems that use more of
its features. Each one shows the same system built in two ways:

- **Command line:** an input file, run with the `packmol` executable, as with Fortran
  Packmol (see [Home](../index.md) for how to install it).
- **Julia:** the same system defined with Julia code, with
  [`structure_type`](@ref), [`PackmolSystem`](@ref) and [`packmol`](@ref) (see
  [Defining systems in Julia](../julia_api/index.md)).

| Tutorial | System |
|---|---|
| [Water and urea](mixture.md) | A mixture of two molecules in a box |
| [Water–chloroform interface](interface.md) | Two solvents in adjacent boxes and a fixed molecule at the interface |
| [Lipid bilayer](bilayer.md) | Oriented molecules, with constraints on some of their atoms |
| [Spherical vesicle](spherical.md) | Lipid layers and water in spherical shells |
| [Solvated protein](solvprotein.md) | A fixed protein in a sphere of water and ions |
| [Protein in glycerol and water](protein_glycerol.md) | A fixed protein in a mixture, in a periodic box |
| [Water drop on graphite](water_drop.md) | A fixed surface, and a drop shaped by constraints |
| [Gases in a zeolite](zeolite.md) | A fixed crystal in its triclinic periodic cell, with gases in its pores |
| [Zika virion](zika.md) | A large system: 2.8 million particles in membranes and spherical shells |

## Getting the files

The input and structure files of the classic examples can be downloaded at once from the Packmol site
([examples.tar.gz](https://m3g.github.io/packmol/examples/examples.tar.gz) or
[examples.zip](https://m3g.github.io/packmol/examples/examples.zip)). Each tutorial links
to the files it needs. Run the examples from the directory that contains the
files: the paths in an input file are relative to the directory of the input file, and
those in a Julia script to the current directory.

The input files of the Packmol site also come in commented versions (`*-comment.inp`),
which are the ones shown in the tutorials.
