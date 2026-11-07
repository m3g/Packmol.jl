# Installation

Users of the legacy (Fortran) Packmol: see [Migrating from legacy Packmol](@ref) for the
differences.

Install the latest Julia version in your system, with [juliaup](https://github.com/JuliaLang/juliaup#juliaup---julia-version-manager)

## Interactive use

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

## Command line interface

Install the `packmol` app with:
```
julia -e 'import Pkg; Pkg.Apps.add("Packmol")'
```

Add the directory of julia binaries (`$HOME/.julia/bin`) to your path, and use packmol with:
```
packmol input.inp
```
(or `packmol -i input.inp`, or `packmol < input.inp`, as Fortran Packmol) as a standalone
application. `-o output.pdb` overrides the output file of the input, `-s N` (or `--seed N`)
the seed of the random number generator, and `-t N` limits the
number of threads used (default: `JULIA_NUM_THREADS` if set, otherwise all the threads;
the input keyword `nthreads` does the same). `--guess-volume` and `--free-volume` estimate
the volume of the system and of the regions of its molecules, without packing it (see
[Volume estimates](@ref)). Without an input file
(and no input redirected), a file browser is opened to choose it. The exit codes are those
of Fortran Packmol: `0` for success, `173` if the packing did not converge (the best
solution is written to the output file, and a solution satisfying only the constraints to
`<output>_FORCED`), `171` for errors in the input, `172` if a file could not be opened,
`174` for command-line errors, and `170` for other errors.

The report of the run goes to the standard output, and can be redirected to a file, as
with Fortran Packmol: `packmol input.inp > packmol.log`. The progress bars, as those of
Fortran Packmol, are lines of `*` characters appended as the packing advances, so the
progress of a long run can be followed in the log (`tail -f packmol.log`). Warnings and
errors are printed to the standard error, and scripts can rely on the exit code. The
global keyword `screen_log no` suppresses all the output of a run (including the
warnings), except the errors.

### Replicas of a system

`-o` writes the packed system to another file than the `output` of the input file, so the
same input file can build several replicas of a system, for independent simulations. Each
replica needs a different seed of the random number generator (with the same seed, every
run gives the same configuration), given with `-s` (or `--seed`). For example, in a Linux
or macOS shell:

```bash
for i in 1 2 3 4 5; do
    packmol -i system.inp -o system_$i.pdb -s $i > system_$i.log
done
```

Each replica can then be built again with its seed. Alternatively, `seed -1` in the input
file (or `-s -1`) takes a new random seed at each run.

The same from Julia:

```julia
using Packmol
for i in 1:5
    system = read_packmol_input("system.inp")
    system.output_file = "system_$i.pdb"
    packmol(system; seed = i)
end
```

### File browser

The file browser can be disabled by setting the environment variable
`PACKMOL_GUI="false"`. It uses tools that come with the operating system (`osascript` on
macOS, PowerShell on Windows, and `zenity` or `kdialog` on Linux); if none of them is
available, the input file is asked for in the terminal.

## Updating

To keep Packmol up-to-date, use:
```
julia> import Pkg; Pkg.update()
```

(or type `] up`, at the `julia>` prompt), and, for the `packmol` app:
```
julia -e 'import Pkg; Pkg.Apps.update("Packmol")'
```
