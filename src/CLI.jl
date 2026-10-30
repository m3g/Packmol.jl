module CLI

import ..Packmol
using LinearAlgebra: BLAS

const usage = """
    Usage: packmol input.inp
           packmol -i input.inp [-o output.pdb]
           packmol < input.inp

    With no input file and no input redirected to stdin, a file dialog is
    opened to choose the input file. Paths in the input file are relative to
    its directory (or to the current directory, if it is read from stdin).

    Options:
      -i <file>      input file
      -o <file>      output file (overrides the input file's `output`)
      -t <n>         maximum number of threads (overrides the input file's
                     `nthreads`; default: JULIA_NUM_THREADS if set, else all)
      -h, --help     show this help
      -v, --version  show the version
    """

# Exit codes of Fortran Packmol (src/exit_codes.f90), kept for compatibility
# with scripts that check them.
const exit_code_success = 0
const exit_code_general_error = 170
const exit_code_input_error = 171
const exit_code_open_file = 172
const exit_code_failed_to_converge = 173
const exit_code_command_line_error = 174

# Whether stdin is a terminal (the libc `isatty`: `stdin` itself is not
# initialized in a trimmed executable).
_stdin_isatty() = Sys.iswindows() ?
    ccall((:_isatty, "msvcrt"), Cint, (Cint,), 0) != 0 :
    ccall(:isatty, Cint, (Cint,), 0) != 0

# Contents of stdin, read with the libc `read` (works for both redirected
# files and pipes, which can't be seeked).
function _read_stdin()
    io = IOBuffer()
    buf = Vector{UInt8}(undef, 65536)
    while true
        n = Sys.iswindows() ?
            Int(ccall((:_read, "msvcrt"), Cint, (Cint, Ptr{UInt8}, Cuint), 0, buf, length(buf))) :
            Int(ccall(:read, Cssize_t, (Cint, Ptr{UInt8}, Csize_t), 0, buf, length(buf)))
        n <= 0 && break
        write(io, view(buf, 1:n))
    end
    return seekstart(io)
end

# The executable is built with all threads (see juliac/build.jl), so JULIA_NUM_THREADS
# does not set their number: it limits the number of threads used, as `-t`. Its
# "auto" and "n,m" (interactive threads) forms are accepted.
function _env_nthreads()
    s = something(Base.get(ENV, "JULIA_NUM_THREADS", nothing), "")
    n = tryparse(Int, first(split(s, ',')))
    return isnothing(n) || n < 1 ? nothing : n
end

function _print_error(@nospecialize(e))
    # The messages are narrowed to `String` (fields typed `AbstractString`)
    msg = e isa ArgumentError ? e.msg :
          e isa ErrorException ? e.msg :
          e isa SystemError ? e.prefix : "unexpected error."
    msg = msg isa String ? strip(msg) : "unexpected error."
    if e isa SystemError
        Packmol._println("ERROR: ", msg, ": ", Libc.strerror(e.errnum))
    else
        Packmol._println("ERROR: ", msg)
    end
end

function _error_exit(@nospecialize(e), code::Int)
    _print_error(e)
    return e isa SystemError ? exit_code_open_file : code
end

function (@main)(ARGS)
    # Print to the `Core` streams (see `src/io/terminal.jl`)
    Packmol._use_core_streams[] = true
    # Packmol makes no BLAS calls: don't let OpenBLAS start its own threads
    BLAS.set_num_threads(1)
    input_file = ""
    output_file = ""
    nthreads = ""
    positional = 0
    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "-h" || arg == "--help"
            Packmol._print(usage)
            return exit_code_success
        elseif arg == "-v" || arg == "--version"
            Packmol._println("packmol ", Packmol._packmol_version)
            return exit_code_success
        elseif (arg == "-i" || arg == "-o" || arg == "-t" || arg == "--threads") && i < length(ARGS)
            if arg == "-i"
                input_file = ARGS[i+1]
            elseif arg == "-o"
                output_file = ARGS[i+1]
            else
                nthreads = ARGS[i+1]
            end
            i += 2
            continue
        elseif !startswith(arg, "-") && positional == 0
            input_file = arg
            positional += 1
        else
            Packmol._println("ERROR: unrecognized command-line argument: ", arg)
            Packmol._print(usage)
            return exit_code_command_line_error
        end
        i += 1
    end
    ntasks = isempty(nthreads) ? 0 : something(tryparse(Int, nthreads), 0)
    if !isempty(nthreads) && ntasks < 1
        Packmol._println("ERROR: the number of threads (-t) must be a positive integer, got: ", nthreads)
        return exit_code_command_line_error
    end
    # Dimension and number type given as types, so that the parser is type-stable
    system = try
        if !isempty(input_file)
            Packmol.read_packmol_input(input_file, Val(3), Float64)
        elseif !_stdin_isatty()
            Packmol.read_packmol_input(_read_stdin(), Val(3), Float64)
        else
            Packmol.read_packmol_input(Packmol.pick_file(), Val(3), Float64)
        end
    catch e
        return _error_exit(e, exit_code_input_error)
    end
    isempty(output_file) || (system.output_file = output_file)
    if ntasks > 0
        system.nthreads = ntasks
    elseif isnothing(system.nthreads)
        system.nthreads = _env_nthreads()
    end
    try
        Packmol.packmol(system)
    catch e
        return _error_exit(e, exit_code_general_error)
    end
    return system.status == :failed_packing ? exit_code_failed_to_converge : exit_code_success
end

end
