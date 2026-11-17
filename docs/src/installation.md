# Installation

Packmol can be installed in these ways:

- **[pip or uv](@ref pip)** and **[binaries](@ref binaries)**: the `packmol` command, a
  standalone executable, with no Julia needed. It runs input files, with every feature
  of them, and is the simplest to integrate in scripts and workflows.
- **[Julia package](@ref julia_package)**: recommended for building systems with an AI
  agent, and to use all the features of Packmol.jl, some of which are only available
  from Julia code.
- **[Julia app](@ref julia_app)**: the `packmol` command, installed with Julia.

All the options run the same packing engine.

## [pip or uv](@id pip)

The `packmol` command, without Julia:

```bash
pip install packmol            # or: uv tool install packmol
packmol input.inp
```

or, without installing it:

```bash
uvx packmol input.inp
```

The `packmol` package of PyPI installs the executable compiled from Packmol.jl from
version 100.0.0 on. Wheels are provided for Linux (x86_64 and aarch64, glibc ≥ 2.28),
macOS (arm64 and x86_64, macOS ≥ 11) and Windows (x86_64). See
[Running Packmol](@ref) for its options.

## [Binaries](@id binaries)

Archives of the executable for each platform, from the
[latest release](https://github.com/m3g/Packmol.jl/releases/latest) (from version 100.0.0
on), or from the [development build](https://github.com/m3g/Packmol.jl/releases/tag/dev)
(the latest commit, not a release):

| System | Architecture | Latest release | Development build |
|:-------|:-------------|:---------------|:------------------|
| Linux | x86_64 | [download](https://github.com/m3g/Packmol.jl/releases/latest/download/packmol-linux-x86_64.tar.gz) | [download](https://github.com/m3g/Packmol.jl/releases/download/dev/packmol-linux-x86_64.tar.gz) |
| Linux | aarch64 (arm64) | [download](https://github.com/m3g/Packmol.jl/releases/latest/download/packmol-linux-aarch64.tar.gz) | [download](https://github.com/m3g/Packmol.jl/releases/download/dev/packmol-linux-aarch64.tar.gz) |
| macOS | Apple silicon (arm64) | [download](https://github.com/m3g/Packmol.jl/releases/latest/download/packmol-macos-arm64.tar.gz) | [download](https://github.com/m3g/Packmol.jl/releases/download/dev/packmol-macos-arm64.tar.gz) |
| macOS | Intel (x86_64) | [download](https://github.com/m3g/Packmol.jl/releases/latest/download/packmol-macos-x86_64.tar.gz) | [download](https://github.com/m3g/Packmol.jl/releases/download/dev/packmol-macos-x86_64.tar.gz) |
| Windows | x86_64 | [download](https://github.com/m3g/Packmol.jl/releases/latest/download/packmol-windows-x86_64.zip) | [download](https://github.com/m3g/Packmol.jl/releases/download/dev/packmol-windows-x86_64.zip) |

Extract the archive anywhere, and run `bin/packmol` (`bin\packmol.exe` on Windows) from
the extracted directory, or add its `bin` directory to the path. Keep the extracted
directory whole: the executable needs the libraries of its `lib` directory. For example,
on Linux:

```bash
curl -LO https://github.com/m3g/Packmol.jl/releases/latest/download/\
packmol-linux-x86_64.tar.gz
tar xzf packmol-linux-x86_64.tar.gz
./packmol-*-linux-x86_64/bin/packmol input.inp
```

The Linux executables need glibc ≥ 2.28 (RHEL/Rocky/Alma 8, Ubuntu 20.04, Debian 11 and
newer). On macOS, archives downloaded with a browser are quarantined: remove the
quarantine with `xattr -dr com.apple.quarantine <extracted directory>`. Each release
lists the SHA256 sums of its archives in `SHA256SUMS.txt`.

See [Running Packmol](@ref) for its options.

## [Julia package](@id julia_package)

The Julia package is the recommended installation for building systems with an **AI
agent**, and the one that supports all the features of Packmol.jl. Besides running input
files, it has:

- the definition of systems in Julia code, without input files (see
  [Defining systems in Julia](julia_api/index.md));
- the **recipes**: solvated systems and membranes built from target densities and
  concentrations (see [Recipes](recipes/index.md));
- user-defined constraints ([`CustomConstraint`](@ref));
- the packed atoms in memory, for analysis with
  [PDBTools.jl](https://github.com/m3g/PDBTools.jl) and the other packages of the Julia
  ecosystem.

An agent that writes a Julia script can, in a single script, obtain and prepare the
molecules, compute the numbers of molecules and the size of the box, pack the system,
and check the result; the script is a record of how the system was built. See
[Using your AI agent](tutorials/ai_agent.md).

Install the latest Julia version in your system, with
[juliaup](https://github.com/JuliaLang/juliaup#juliaup---julia-version-manager). Then
install the Packmol package within Julia:

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

To keep Packmol up-to-date, use:
```julia-repl
julia> import Pkg; Pkg.update()
```
(or type `] up`, at the `julia>` prompt).

## [Julia app](@id julia_app)

The `packmol` command can also be installed with Julia, on every system where Julia
runs, including those for which there are no binaries (e.g. Windows on ARM):

```bash
julia -e 'import Pkg; Pkg.Apps.add("Packmol")'
```

Add the directory of the Julia binaries (`$HOME/.julia/bin`) to your path. To update it:

```bash
julia -e 'import Pkg; Pkg.Apps.update("Packmol")'
```

## Legacy Packmol

Packmol.jl reads the input files of the legacy (Fortran) Packmol, and its command line is
the same. See [Migrating from legacy Packmol](@ref) for the differences. The legacy
Packmol remains available at [github.com/m3g/packmol](https://github.com/m3g/packmol),
and from PyPI, with `pip install "packmol<100"` (the versions before 100).
