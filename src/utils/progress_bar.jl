#
# A progress bar printed as Fortran Packmol's: a header line, then the bar, drawn by
# appending `*` characters to a line of `stdout`, ended by the elapsed time:
#
#       Iterations: |0                                                     100%|
#                   |************************************************************| 0.3 s
#
# Nothing is redrawn (no carriage returns or escape sequences), so the bar is the same on
# a terminal, in a file (`packmol input.inp > log`), in a pipe and in a notebook, and the
# progress of a long loop can be followed in a log (`tail -f`). (ProgressMeter.jl, which
# redraws in place, does not compile with `juliac --trim`.) `next!` is thread-safe.
#
mutable struct ProgressBar
    const n::Int
    const indent::Int # length of the label, before the bar
    const width::Int
    const tstart::Float64
    const lock::ReentrantLock
    counter::Int
    nprinted::Int # characters of the bar already printed
    finished::Bool
end
function ProgressBar(n::Integer; desc::String="Progress: ", width::Int=60)
    p = ProgressBar(n, length(desc), width, time(), ReentrantLock(), 0, 0, false)
    _with_stdout() do io
        # The header: `0` at the start of the bar, `100%` at its end
        print(io, desc, "|0")
        _print_repeated(io, ' ', width - 5)
        print(io, "100%|\n")
        _print_repeated(io, ' ', p.indent)
        print(io, '|')
        flush(io)
    end
    return p
end

_format_seconds(t::Real) = t < 60 ? @sprintf("%.1f s", t) :
    @sprintf("%d:%02d:%02d", round(Int, t) ÷ 3600, (round(Int, t) % 3600) ÷ 60, round(Int, t) % 60)

# Append the characters of the bar up to the fraction `p.counter / p.n` of its width
function _extend!(p::ProgressBar)
    target = p.n > 0 ? min(p.width, (p.width * p.counter) ÷ p.n) : p.width
    if target > p.nprinted
        _with_stdout() do io
            _print_repeated(io, '*', target - p.nprinted)
            flush(io)
        end
        p.nprinted = target
    end
    return p
end

function next!(p::ProgressBar)
    @lock p.lock begin
        p.finished && return p
        p.counter += 1
        _extend!(p)
    end
    return p
end

# Completes the bar (also when stopped early, e.g. on convergence), with the elapsed time
function finish!(p::ProgressBar)
    @lock p.lock begin
        p.finished && return p
        p.counter = p.n
        _extend!(p)
        p.finished = true
        _println("| ", _format_seconds(time() - p.tstart))
    end
    return p
end
