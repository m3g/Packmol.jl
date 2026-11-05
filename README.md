[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://m3g.github.io/Packmol.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://m3g.github.io/Packmol.jl/dev)
[![Tests](https://img.shields.io/badge/build-passing-green)](https://github.com/m3g/Packmol.jl/actions)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# Packmol.jl

Packmol creates initial configurations for molecular dynamics simulations, by packing
molecules in regions of space defined by geometric constraints, such that the atoms of
different molecules are kept at a safe distance from each other.

Packmol.jl is the new implementation of [Packmol](http://github.com/m3g/packmol). It reads
the input files of the legacy (Fortran) Packmol and runs them with the same command line,
and it can also be used from Julia code, where systems are defined without input files and
can be built from target densities and concentrations (recipes). See the
[documentation](https://m3g.github.io/Packmol.jl/stable), and
[Migrating from legacy Packmol](https://m3g.github.io/Packmol.jl/stable/migrating/) for the
differences.

## Installation

Install the latest Julia version in your system, with
[juliaup](https://github.com/JuliaLang/juliaup#juliaup---julia-version-manager).

### Command line

Install the `packmol` app with:
```
julia -e 'import Pkg; Pkg.Apps.add("Packmol")'
```

Add the directory of julia binaries (`$HOME/.julia/bin`) to your path, and run:
```
packmol input.inp
```
or `packmol < input.inp`, as the legacy Packmol. `-o output.pdb` overrides the output file
of the input, `-s N` its random seed, and `-t N` limits the number of threads; `packmol -h`
lists the options. Without an input file (and no input redirected), a file browser is
opened to choose it. The exit codes are those of the legacy Packmol: `0` for success,
`173` if the packing did not converge (the best solution is written to the output file,
and a solution satisfying only the constraints to `<output>_FORCED`), `171` for errors in
the input, `172` if a file could not be opened, `174` for command-line errors, and `170`
for other errors.

The file browser uses tools that come with the operating system (`osascript` on macOS,
PowerShell on Windows, and `zenity` or `kdialog` on Linux), and it can be disabled by
setting the environment variable `PACKMOL_GUI="false"`.

### Julia package

```julia-repl
julia> import Pkg; Pkg.add("Packmol")

julia> using Packmol

julia> packmol("input.inp")
```

`packmol()`, without an input file, opens the file browser. Systems can also be defined
in Julia, without input files, and built from densities and concentrations with the
recipes: see the documentation.

### Updating

```
julia -e 'import Pkg; Pkg.Apps.update("Packmol"); Pkg.update()'
```

### Legacy (Fortran) Packmol

The legacy Packmol is still available from [m3g/packmol](http://github.com/m3g/packmol),
and it can be run from Julia, with `run_packmol_legacy`, which works as `packmol`. It
requires the `Packmol_jll` package, which provides the Fortran binary:

```julia-repl
julia> import Pkg; Pkg.add("Packmol_jll")

julia> using Packmol, Packmol_jll

julia> run_packmol_legacy()
```

## References

If Packmol contributed to your work, please cite the following references:

- L. Martínez, R. Andrade, E. G. Birgin, J. M. Martínez. Packmol: A package for building
  initial configurations for molecular dynamics simulations. *Journal of Computational
  Chemistry*, 30(13):2157-2164, 2009. [doi:10.1002/jcc.21224](https://doi.org/10.1002/jcc.21224)

- J. M. Martínez and L. Martínez. Packing optimization for automated generation of complex
  system's initial configurations for molecular dynamics and docking. *Journal of
  Computational Chemistry*, 24(7):819-825, 2003. [doi:10.1002/jcc.10216](https://doi.org/10.1002/jcc.10216)
