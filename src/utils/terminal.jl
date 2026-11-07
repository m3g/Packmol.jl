#
# Terminal input and output that compiles with `juliac --trim`.
#
# `stdout`, `stderr` and `stdin` are non-constant globals of abstract type `IO`, so
# calls that print to (or read from) them cannot be resolved at compile time, as
# `--trim` requires. These helpers branch on the concrete types the streams usually
# have (a terminal, a pipe, a file), falling back to the `Core` streams, which have
# concrete types, otherwise (e.g. the output of an IJulia notebook cell).
#
# In an executable compiled with `juliac --trim`, `stdout` and `stderr` are not
# initialized at all (they keep the streams of the compilation), so the command line
# app (`CLI.jl`) sets `_use_core_streams[] = true`, to print to the `Core` streams only.
#
const _use_core_streams = Ref(false)

# (A macro, rather than a function taking the stream, because passing the abstractly
# typed stream as an argument would already be a call that cannot be resolved.)
macro _with_stream(f, stream, fallback)
    quote
        io = $(esc(stream))
        if io isa Base.TTY
            $(esc(f))(io)
        elseif io isa Base.PipeEndpoint
            $(esc(f))(io)
        elseif io isa IOStream
            $(esc(f))(io)
        elseif io isa Base.DevNull
            $(esc(f))(io)
        else
            $(esc(f))($(esc(fallback)))
        end
    end
end
function _with_stdout(f::F) where {F}
    _use_core_streams[] && return f(Core.stdout)
    return @_with_stream(f, stdout, Core.stdout)
end
function _with_stderr(f::F) where {F}
    _use_core_streams[] && return f(Core.stderr)
    return @_with_stream(f, stderr, Core.stderr)
end

# Whether `stderr` is a terminal (for the progress bar to be redrawn in place).
function _stderr_is_terminal()
    _use_core_streams[] || return stderr isa Base.TTY
    return @static if Sys.iswindows()
        ccall(:_isatty, Cint, (Cint,), 2) != 0
    else
        ccall(:isatty, Cint, (Cint,), 2) != 0
    end
end

# Print the elements of the tuple `xs` one by one, so that each `print` call has a
# concrete argument type (printing a heterogeneous tuple at once does not). `map`
# over a tuple is unrolled, unlike a loop over it, or a recursion, which inference
# eventually widens. (A tuple argument, rather than varargs, because a varargs
# method called with splatted arguments is not always specialized on them.)
_print_each(io::IO, xs::Tuple) = (map(x -> print(io, x), xs); nothing)

_print(xs...) = _with_stdout(io -> _print_each(io, xs))
_println(xs...) = _with_stdout(io -> _print_each(io, (xs..., '\n')))
_flush_stdout() = _with_stdout(flush)

# Print `x` `n` times (`repeat`, `lpad` and `^` on strings do not compile with `--trim`).
function _print_repeated(io::IO, x, n::Integer)
    for _ in 1:n
        print(io, x)
    end
end

# `@printf` to `stdout`, through `_print`.
macro _printf(format, args...)
    # (The format must reach `@sprintf` as a literal string, so it is not escaped.)
    sprintf = Expr(:macrocall, GlobalRef(@__MODULE__, Symbol("@sprintf")), __source__, format, map(esc, args)...)
    return :($(GlobalRef(@__MODULE__, :_print))($sprintf))
end

# Read a line from the standard input, with the C library's `read`, because reading
# from Julia's `stdin` (a libuv stream) does not compile with `--trim`.
function _readline()
    _flush_stdout()
    line = IOBuffer()
    byte = Ref{UInt8}(0)
    while true
        n = @static if Sys.iswindows()
            ccall(:_read, Cint, (Cint, Ptr{UInt8}, Cuint), 0, byte, 1)
        else
            ccall(:read, Cssize_t, (Cint, Ptr{UInt8}, Csize_t), 0, byte, 1)
        end
        (n <= 0 || byte[] == UInt8('\n')) && break
        write(line, byte[])
    end
    return String(take!(line))
end
