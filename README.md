# Packmol

Currently this serves as a multi-platform runner for Packmol.

## Installation

1. Install Julia in your system, with [juliaup](https://github.com/JuliaLang/juliaup#juliaup---julia-version-manager)

2. Launch julia and install this package:

```julia
julia> import Pkg; Pkg.add("Packmol")
```

## Usage

### Interactive use

Run packmol in your input file with:

```julia
julia> using Packmol

julia> run_packmol()
```

This will open a file browser, from which you can choose the input file for `packmol`. 
Packmol will run immediately once the file is open.

Alternativelly, you can provide the path to the file explicitly, with:

```julia
julia> run_packmol(raw"C:\users\my_user\my_files\my_input_file.inp")
```

### Scripting

For command-line usage, we suggest the following procedure:

Start Julia and create an environment where packmol is installed:

```julia
julia> import Pkg

julia> Pkg.activate("Packmol", shared=true)

julia> Pkg.add("Packmol")

julia> exit()
```

Next, write a script like:

```julia
import Pkg; Pkg.activate("Packmol", shared=true) # Activate Packmol environment
using Packmol
input_file = raw"C:\users\my_user\my_files\my_input_file.inp" 
run_packmol(input_file)
```

This script can be executed from the command line with
```
julia script.jl
```

If you want to run `packmol` in multiple files, a vector of files
can be provided:

```julia
import Pkg; Pkg.activate("Packmol", shared=true) # Activate Packmol environment
using Packmol
input_files = [
    raw"C:\users\my_user\my_files\my_input_file1.inp",
    raw"C:\users\my_user\my_files\my_input_file2.inp", 
    raw"C:\users\my_user\my_files\my_input_file3.inp",
]
for input_file in input_files
    run_packmol(input_file)
end
```




