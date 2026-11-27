# Running Packmol

The `packmol` command (installed with pip or uv, from the binaries, or as the Julia app;
see [Installation](installation.md)) runs an input file:

```bash
packmol input.inp [-o output.pdb]
```

The paths in the input file (structures, `output`, `restart_from`) are relative to its
directory, so it runs from any directory.

Options:

- `-i input.inp`: the input file (the same as `packmol input.inp`);
- `-o output.pdb`: overrides the output file of the input;
- `-s N` (or `--seed N`): the seed of the random number generator;
- `-t N`: the number of threads used (default: `JULIA_NUM_THREADS` if set, otherwise all
  the threads; the input keyword `nthreads` does the same);
- `--guess-volume`, `--free-volume` and `--molecular-volume`: estimate the volume of
  the system, of the regions of its molecules, and of the molecules themselves, without
  packing it (see [Volume estimates](@ref)).

Without an input file, a file browser is opened to choose it.

The exit codes are `0` for success, `173` if the packing did not converge (the best
solution is written to the output file, and a solution satisfying only the constraints to
`<name>-FORCED.<ext>`, e.g. `system-FORCED.pdb` for `system.pdb`), `171` for errors in the
input, `172` if a file could not be opened, `174` for command-line errors, and `170` for
other errors.

The report of the run goes to the standard output, and can be redirected to a file:
`packmol input.inp > packmol.log`. The progress bars are lines of `*` characters appended
as the packing advances, so the progress of a long run can be followed in the log
(`tail -f packmol.log`). Warnings and errors are printed to the standard error, and
scripts can rely on the exit code. The global keyword `screen_log no` suppresses all the
output of a run (including the warnings), except the errors.

## Replicas of a system

`-o` writes the packed system to another file than the `output` of the input file, so the
same input file can build several replicas of a system, for independent simulations. Each
replica needs a different seed of the random number generator (with the same seed, every
run gives the same configuration), given with `-s` (or `--seed`). For example, in a Linux
or macOS shell:

```bash
for i in 1 2 3 4 5; do
    packmol system.inp -o system_$i.pdb -s $i > system_$i.log
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

## File browser

The file browser can be disabled by setting the environment variable
`PACKMOL_GUI="false"`. It uses tools that come with the operating system (`osascript` on
macOS, PowerShell on Windows, and `zenity` or `kdialog` on Linux); if none of them is
available, the input file is asked for in the terminal.

## Input redirected with `<`

`packmol < input.inp`, the way the legacy Packmol was run, still works, for
compatibility. Prefer `packmol input.inp`: with the input redirected, `packmol` does not
know where the input file is, so the paths in it (structures, `output`) are relative to
the current directory: it must be run from the directory of the input file.
