[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://m3g.github.io/Packmol.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://m3g.github.io/Packmol.jl/dev)
[![Tests](https://img.shields.io/badge/build-passing-green)](https://github.com/m3g/Packmol.jl/actions)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# Packmol.jl

This Julia package is a new implementation of [packmol](http://github.com/m3g/packmol),
which can also run the classical Fortran program (see [Legacy (Fortran) Packmol](#legacy-fortran-packmol)).

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
(or `packmol -i input.inp`) as a standalone application. Without an input file, a file
browser is opened to choose it.

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






