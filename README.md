[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://m3g.github.io/Packmol.jl/stable)
[![Tests](https://img.shields.io/badge/build-passing-green)](https://github.com/m3g/Packmol.jl/actions)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# Packmol.jl

Currently this Julia package serves as a multi-platform runner for [packmol](http://github.com/m3g/packmol). In the future this will be an 
independent package with a faster and improved version of the Packmol package.

## Installation

1. Install Julia in your system, with [juliaup](https://github.com/JuliaLang/juliaup#juliaup---julia-version-manager)

2. Launch julia and install this package:

```julia-repl
julia> import Pkg; Pkg.add("Packmol")
```

## Usage

### Interactive use

Run packmol in your input file with:

```julia-repl
julia> using Packmol

julia> run_packmol()
```

This will open a file browser, from which you can choose the input file for `packmol`. 
Packmol will run immediately once the file is open.

Alternativelly, you can provide the path to the file explicitly, with:

```julia-repl
julia> run_packmol(raw"C:\users\my_user\my_files\my_input_file.inp")
```

### Scripting

For command-line usage, we suggest the following procedure:

Start Julia and create an environment where packmol is installed:

```julia-repl
julia> import Pkg

julia> Pkg.activate("Packmol", shared=true)

julia> Pkg.add("Packmol")

julia> exit()
```

Next, copy this into a `packmol.jl` file:

```julia
import Pkg; Pkg.activate("Packmol", shared=true) # Activate Packmol environment
using Packmol
for input_file in ARGS
    run_packmol(input_file)
end
```

This script can be executed from the command line with
```
julia packmol.jl packmol_input1.inp packmol_input2.inp
```

to run `packmol` for one or more input files.

To keep Packmol up-to-date, use:
```julia-repl
julia> import Pkg; Pkg.activate("Packmol"; shared=true); Pkg.update()
```

(or type `] activate @Packmol` and then `] up`, at the `julia>` prompt). 




