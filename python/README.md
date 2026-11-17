# Packmol

Packmol creates initial configurations for molecular dynamics simulations, by packing
molecules in defined regions of space. The packing guarantees that short-range repulsive
interactions do not disrupt the simulations.

This package installs the `packmol` command:

```bash
pip install packmol      # or: uv tool install packmol
packmol input.inp [-o output.pdb]
uvx packmol input.inp    # without installing it
```

Since version 100.0.0, the executable is compiled from
[Packmol.jl](https://github.com/m3g/Packmol.jl), the Julia implementation of Packmol, which
reads the same input files. See the
[documentation](https://m3g.github.io/Packmol.jl), and its page on migrating from the
Fortran version. The Fortran Packmol (versions 21.x) remains available at
[github.com/m3g/packmol](https://github.com/m3g/packmol) (`pip install "packmol<100"`).

Please cite: L. Martínez, R. Andrade, E. G. Birgin, J. M. Martínez, *PACKMOL: A package for
building initial configurations for molecular dynamics simulations*, J. Comput. Chem.
30, 2157-2164 (2009), doi:[10.1002/jcc.21224](https://doi.org/10.1002/jcc.21224).
