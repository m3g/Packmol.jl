#
# A minimal progress bar, replacing ProgressMeter.jl, which does not compile with
# `juliac --trim`. Printed to `stderr`, as ProgressMeter's, and redrawn at most
# every `dt` seconds. On a terminal the bar (and the `showvalues` lines below it)
# is redrawn in place; otherwise it is printed only once, by `finish!`, so that
# logs are not flooded. `next!` is thread-safe.
#
mutable struct ProgressBar
    const n::Int
    const desc::String
    const barlen::Int
    const dt::Float64
    const tstart::Float64
    const lock::ReentrantLock
    counter::Int
    tlast::Float64
    showvalues::String # the last `showvalues`, printed below the bar
    nlines::Int # lines printed by the last redraw, to be overwritten by the next one
    finished::Bool
end
ProgressBar(n::Integer; desc::String="Progress: ", barlen::Int=40, dt::Real=0.1) =
    ProgressBar(n, desc, barlen, dt, time(), ReentrantLock(), 0, 0.0, "", 0, false)

function _format_seconds(t::Real)
    s = round(Int, t)
    return @sprintf("%d:%02d:%02d", s ÷ 3600, (s % 3600) ÷ 60, s % 60)
end

# `showvalues` is a tuple of `(label, value)` tuples, printed one per line.
function _set_showvalues!(p::ProgressBar, showvalues::Tuple)
    isempty(showvalues) && return p
    buf = IOBuffer()
    # (`map` rather than a loop, for each value to be printed with its own type.)
    map(((label, value),) -> _print_each(buf, ("\n  ", label, ":  ", value)), showvalues)
    p.showvalues = String(take!(buf))
    return p
end

function _draw(p::ProgressBar)
    frac = p.n > 0 ? min(p.counter / p.n, 1.0) : 1.0
    nfilled = round(Int, frac * p.barlen)
    elapsed = time() - p.tstart
    buf = IOBuffer()
    print(buf, p.desc, @sprintf("%3d", round(Int, 100 * frac)), "%|")
    _print_repeated(buf, '█', nfilled)
    _print_repeated(buf, ' ', p.barlen - nfilled)
    print(buf, "| ")
    if p.finished
        print(buf, "Time: ", _format_seconds(elapsed))
    elseif p.counter > 0
        print(buf, "ETA: ", _format_seconds(elapsed * (p.n - p.counter) / p.counter))
    end
    print(buf, p.showvalues)
    str = String(take!(buf))
    terminal = _stderr_is_terminal()
    _with_stderr() do io
        # Move back to the beginning of the previous redraw, clearing it.
        if terminal && p.nlines > 0
            print(io, "\r\e[K")
            _print_repeated(io, "\e[1A\r\e[K", p.nlines - 1)
        end
        print(io, str)
        p.finished && print(io, '\n')
        flush(io)
    end
    p.nlines = 1 + count(==('\n'), p.showvalues)
    return p
end

function next!(p::ProgressBar; showvalues::Tuple=())
    @lock p.lock begin
        p.counter += 1
        _set_showvalues!(p, showvalues)
        t = time()
        if _stderr_is_terminal() && t - p.tlast >= p.dt
            p.tlast = t
            _draw(p)
        end
    end
    return p
end

function finish!(p::ProgressBar; showvalues::Tuple=())
    @lock p.lock begin
        _set_showvalues!(p, showvalues)
        if !p.finished
            # Complete, as ProgressMeter's, even if stopped early (e.g. on convergence)
            p.counter = p.n
            p.finished = true
            _draw(p)
        end
    end
    return p
end
