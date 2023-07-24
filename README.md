# Packmol

Currently this serves as a multi-platform runner for Packmol.

## Installation

1. Install Julia in your system, with [juliaup](https://github.com/JuliaLang/juliaup#juliaup---julia-version-manager)

2. Launch julia and install this package:

```julia
julia> import Pkg; Pkg.add("Packmol")
```

## Usage

3. Run packmol in your input file with:

```julia
julia> using Packmol

julia> run_packmol("C:\\users\\my_user\\my_files\\my_input_file.inp")
```

**Note:** the `\\` are necessary to specify paths on Windows systems, because `\` is a scape-character.
