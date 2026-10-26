export run_packmol_legacy

"""
    run_packmol_legacy()
    run_packmol_legacy(input_file::String)

Runs the packmol executable with the input file `input_file`. This will run
the classical `http://m3g.iqm.unicamp.br/packmol` program, which is a
pre-compiled binary. The input file is a text file with the same syntax as
the packmol input files.

The binary is provided by the `Packmol_jll` package, which must be installed
and loaded (`using Packmol_jll`) for `run_packmol_legacy` to be available.

If no input file is provided, a file explorer will be opened to choose the
input file (see [`Packmol.pick_file`](@ref)).

To disable the file explorer, set the environment variable `PACKMOL_GUI="false"`.

"""
function run_packmol_legacy end

# Fallback for when the `PackmolJLLExt` extension is not loaded: its own
# `run_packmol_legacy()`/`run_packmol_legacy(::String)` methods are more
# specific than this one.
function run_packmol_legacy(args...)
    throw(ArgumentError("""\n
        `run_packmol_legacy` runs the legacy Fortran packmol binary, provided by the Packmol_jll package.
        Install it with `import Pkg; Pkg.add("Packmol_jll")` and load it with `using Packmol_jll`.

    """))
end
