# Packmol.jl

Packmol creates initial configurations for molecular dynamics simulations, by packing
molecules in regions of space defined by geometric constraints, such that the atoms of
different molecules are kept at a safe distance from each other.

Packmol.jl is the Julia implementation of [Packmol](http://github.com/m3g/packmol). It
reads the input files of Fortran Packmol and runs them with the same command line, and it
can also be used from Julia code, where systems are defined without input files and can
be built from target densities and concentrations (see [Recipes](recipes.md)). It uses
[CellListMap.jl](https://github.com/m3g/CellListMap.jl) for the computation of distances
and [SPGBox.jl](https://github.com/m3g/SPGBox.jl) as the optimizer.

## Installation

Install the latest Julia version in your system, with [juliaup](https://github.com/JuliaLang/juliaup#juliaup---julia-version-manager)

### Interactive use

Install the Packmol package within Julia:

```julia-repl
julia> import Pkg; Pkg.add("Packmol")
```

Run packmol in your input file with:

```julia-repl
julia> using Packmol

julia> packmol()
```

This will open a file browser, from which you can choose the input file for `packmol`.
Packmol will run immediately once the file is open.

Alternatively, you can provide the path to the file explicitly, with:

```julia-repl
julia> packmol(raw"C:\users\my_user\my_files\my_input_file.inp")
```

### Command line interface

Install the `packmol` app with:
```
julia -e 'import Pkg; Pkg.Apps.add("Packmol")'
```

Add the directory of julia binaries (`$HOME/.julia/bin`) to your path, and use packmol with:
```
packmol input.inp
```
(or `packmol -i input.inp`, or `packmol < input.inp`, as Fortran Packmol) as a standalone
application. `-o output.pdb` overrides the output file of the input, and `-t N` limits the
number of threads used (default: `JULIA_NUM_THREADS` if set, otherwise all the threads;
the input keyword `nthreads` does the same). Without an input file
(and no input redirected), a file browser is opened to choose it. The exit codes are those
of Fortran Packmol: `0` for success, `173` if the packing did not converge (the best
solution is written to the output file, and a solution satisfying only the constraints to
`<output>_FORCED`), `171` for errors in the input, `172` if a file could not be opened,
`174` for command-line errors, and `170` for other errors.

The file browser can be disabled by setting the environment variable
`PACKMOL_GUI="false"`. It uses tools that come with the operating system (`osascript` on
macOS, PowerShell on Windows, and `zenity` or `kdialog` on Linux); if none of them is
available, the input file is asked for in the terminal.

### Updating

To keep Packmol up-to-date, use:
```
julia> import Pkg; Pkg.update()
```

(or type `] up`, at the `julia>` prompt).

## Documentation

- [Tutorials](tutorials/index.md): complete examples, each with an input file for the
  command line and the same system defined in Julia.
- [Input files](input_files.md): the global keywords and the `structure ... end structure`
  blocks.
- [Constraints](constraints.md): the geometric regions (`box`, `sphere`, `plane`,
  `cylinder`, `ellipsoid`, ...) in which molecules, or some of their atoms, are packed.
- [Periodic boundary conditions](periodic_boundary_conditions.md): orthorhombic,
  triclinic, rhombic dodecahedral and truncated octahedral cells.
- [Defining systems in Julia](julia_api.md): building and packing a system with
  `structure_type`, `PackmolSystem` and `packmol`, without an input file.
- [Recipes](recipes.md): solvated systems and membranes built from target densities and
  concentrations, instead of numbers of molecules.
- [Concentration units](concentration_units.md): the conversions between concentration
  units (`cconvert`) used by the recipes.
- [Reference](reference.md): the docstrings of the public functions.

## Legacy (Fortran) Packmol

The classical Fortran program can still be run from Julia, with
[`run_packmol_legacy`](@ref), which works as `packmol`. It requires the `Packmol_jll`
package, which provides the Fortran binary:

```julia-repl
julia> import Pkg; Pkg.add("Packmol_jll")

julia> using Packmol, Packmol_jll

julia> run_packmol_legacy()
```
