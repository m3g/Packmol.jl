# Packmol.jl

Currently this Julia package serves as a multi-platform runner for [packmol](http://github.com/m3g/packmol). In the future this will be an 
independent package with a faster and improved version of the Packmol package.

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

julia> run_packmol()
```

This will open a file browser, from which you can choose the input file for `packmol`. 
Packmol will run immediately once the file is open.

Alternatively, you can provide the path to the file explicitly, with:

```julia-repl
julia> run_packmol(raw"C:\users\my_user\my_files\my_input_file.inp")
```

### Command line interface

Install the `packmol` app with:
```
julia -e 'import Pkg; Pkg.Apps.add("Packmol")'
```

Add the directory of julia binaries (`$HOME/.julia/bin`) to your path, and use packmol with:
```
packmol -i input.inp
```
as a standalone application.

### Updating

To keep Packmol up-to-date, use:
```
julia> import Pkg; Pkg.update()
```

(or type `] up`, at the `julia>` prompt). 

Additionally, it is possible to disable the loading of the file-dialog machinery
by setting the system environment variable `PACKMOL_GUI="false"`. This might be
important to run packmol through this interface in computers without a GUI. 

## See also

[MolSimToolkit.jl/Packmol-Input-Creator](https://m3g.github.io/MolSimToolkit.jl/stable/system_setup/#Packmol-Input-Creator)

# Reference

```@docs
run_packmol
```

```@docs
pack_monoatomic!
```

