[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://m3g.github.io/Packmol.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://m3g.github.io/Packmol.jl/dev)
[![CI](https://github.com/m3g/Packmol.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/m3g/Packmol.jl/actions/workflows/CI.yml)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# Packmol.jl

Packmol creates initial configurations for molecular dynamics simulations, by packing
molecules in regions of space defined by geometric constraints, such that the atoms of
different molecules are kept at a safe distance from each other.

Packmol.jl is the new implementation of [Packmol](https://doi.org/10.1002/jcc.21224). It
reads the input files of the legacy (Fortran) Packmol, is much faster and multithreaded
(systems of millions of atoms pack in minutes), and adds periodic cells of several shapes,
new constraints, and a Julia API with recipes for solvated systems and membranes.

**[Documentation](https://m3g.github.io/Packmol.jl/stable)**: the
[installation](https://m3g.github.io/Packmol.jl/stable/installation/) (pip or uv,
standalone binaries, or the Julia package), [running
Packmol](https://m3g.github.io/Packmol.jl/stable/running/), the tutorials, and
[migrating from the legacy
Packmol](https://m3g.github.io/Packmol.jl/stable/migrating/). The changes of each version
are in the [CHANGELOG](CHANGELOG.md). To contribute, see [CONTRIBUTING](CONTRIBUTING.md).

## Quick start

```bash
pip install packmol            # or: uv tool install packmol
packmol input.inp [-o output.pdb]
```

or, in Julia:

```julia-repl
julia> import Pkg; Pkg.add("Packmol")

julia> using Packmol

julia> packmol("input.inp")
```

## References

If Packmol contributed to your work, please cite the following references:

- L. Martínez, R. Andrade, E. G. Birgin, J. M. Martínez. Packmol: A package for building
  initial configurations for molecular dynamics simulations. *Journal of Computational
  Chemistry*, 30(13):2157-2164, 2009. [doi:10.1002/jcc.21224](https://doi.org/10.1002/jcc.21224)

- J. M. Martínez and L. Martínez. Packing optimization for automated generation of complex
  system's initial configurations for molecular dynamics and docking. *Journal of
  Computational Chemistry*, 24(7):819-825, 2003. [doi:10.1002/jcc.10216](https://doi.org/10.1002/jcc.10216)
