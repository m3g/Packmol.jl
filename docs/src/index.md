# Packmol.jl

This Julia package is a new implementation of [packmol](http://github.com/m3g/packmol),
which can also run the classical Fortran program (see [Legacy (Fortran) Packmol](@ref)).

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

### Legacy (Fortran) Packmol

`packmol` runs the new Julia implementation of Packmol. The classical Fortran
program can still be run, with `run_packmol_legacy`, which works in the same way.
It requires the `Packmol_jll` package, which provides the Fortran binary:

```julia-repl
julia> import Pkg; Pkg.add("Packmol_jll")

julia> using Packmol, Packmol_jll

julia> run_packmol_legacy()
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

### Updating

To keep Packmol up-to-date, use:
```
julia> import Pkg; Pkg.update()
```

(or type `] up`, at the `julia>` prompt). 

Additionally, it is possible to disable the file browser by setting the system
environment variable `PACKMOL_GUI="false"`. The file browser uses tools that come with
the operating system (`osascript` on macOS, PowerShell on Windows, and `zenity` or
`kdialog` on Linux); if none of them is available, the input file is asked for in the
terminal.

## See also

[Recipes](experimental/recipes.md): building solvated systems from target densities and
concentrations (formerly the Packmol Input Creator of MolSimToolkit.jl).

# Reference

See [`packmol`](@ref) for the Julia implementation.

```@docs
run_packmol_legacy
Packmol.pick_file
```

