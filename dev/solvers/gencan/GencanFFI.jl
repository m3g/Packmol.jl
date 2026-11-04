#
# GENCAN (the optimizer of Fortran Packmol) called from Julia, for the solver experiments
# of dev/solvers/ only: never a dependency of Packmol.jl. Restored from the bindings
# removed in e07c222 (src/gencan/), as a standalone module.
#
# The library is compiled with gfortran on first use, in dev/solvers/gencan/build/.
#
module GencanFFI

using Libdl

const FORTRAN_DIR = joinpath(@__DIR__, "fortran")
const SOURCES = joinpath.(FORTRAN_DIR, (
    "callback_bridge.f90", "gencan_ieee_signal_routines.f90", "gencan_core.f", "ahestetic.f90",
))
const LIBPATH = joinpath(@__DIR__, "build", "libgencan." * Libdl.dlext)
const LIB = Ref{Ptr{Cvoid}}(C_NULL)

function libhandle()
    LIB[] != C_NULL && return LIB[]
    if !isfile(LIBPATH) || any(s -> mtime(s) > mtime(LIBPATH), SOURCES)
        mkpath(dirname(LIBPATH))
        modcache = mktempdir()
        run(`gfortran -shared -fPIC -O2 -o $LIBPATH $SOURCES -J $modcache`)
    end
    LIB[] = Libdl.dlopen(LIBPATH)
    return LIB[]
end

# The Fortran callbacks have no user-data slot: the current objective and the
# per-evaluation hook are kept in globals (calls are never reentrant).
const CURRENT_FG = Ref{Any}(nothing)
const CURRENT_ON_EVAL = Ref{Any}(nothing)

# `on_eval(x, g, f)` is called after every evaluation; returning `true` aborts GENCAN
# (the callback reports an error to it), the way to stop it early.
function fg_callback(n::Cint, x_ptr::Ptr{Cdouble}, f_ptr::Ptr{Cdouble}, g_ptr::Ptr{Cdouble})::Cint
    try
        x = unsafe_wrap(Array, x_ptr, Int(n); own=false)
        g = unsafe_wrap(Array, g_ptr, Int(n); own=false)
        f = CURRENT_FG[](g, x)
        unsafe_store!(f_ptr, f)
        on_eval = CURRENT_ON_EVAL[]
        stop = on_eval === nothing ? false : on_eval(x, g, f)
        return Cint(stop ? 1 : 0)
    catch e
        @error "error in the objective called by GENCAN" exception = (e, catch_backtrace())
        return Cint(1)
    end
end

"""
    gencan!(fg!, x; lower, upper, nitmax, nfevalmax, epsgpsn, on_eval) -> NamedTuple

Minimizes `fg!(g, x) -> f` from `x` (updated in place) with GENCAN (`easygencan`), within
the bounds `lower`/`upper` (`nothing`: unbounded). Returns
`(; f, gpsupn, iter, fcnt, gcnt, cgcnt, inform)`; `inform < 0` when `on_eval` stopped it.
"""
function gencan!(
    fg!, x::Vector{Float64};
    lower=nothing, upper=nothing,
    nitmax::Int=20, nfevalmax::Int=1000, epsgpsn::Float64=1e-6, on_eval=nothing,
)
    n = length(x)
    CURRENT_FG[] = fg!
    CURRENT_ON_EVAL[] = on_eval
    try
        h = libhandle()
        cfun = @cfunction(fg_callback, Cint, (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}))
        ccall(Libdl.dlsym(h, :set_fg_callback), Cvoid, (Ptr{Cvoid},), cfun)
        l = isnothing(lower) ? fill(-1.0e20, n) : Vector{Float64}(lower)
        u = isnothing(upper) ? fill(1.0e20, n) : Vector{Float64}(upper)
        # GENCAN does not accept infinite bounds
        l .= max.(l, -1.0e20)
        u .= min.(u, 1.0e20)
        g = zeros(n)
        wi = zeros(Cint, n)
        wd = zeros(8n)
        f = Ref(0.0); gpsupn = Ref(0.0)
        iter = Ref{Cint}(0); fcnt = Ref{Cint}(0); gcnt = Ref{Cint}(0); cgcnt = Ref{Cint}(0)
        inform = Ref{Cint}(0)
        ccall(
            Libdl.dlsym(h, :easygencan_), Cvoid,
            (
                Ref{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},  # n, x, l, u
                Ref{Cint}, Ptr{Cdouble}, Ptr{Cdouble},                # m, lambda, rho
                Ref{Cdouble}, Ref{Cint}, Ref{Cint},                   # epsgpsn, maxit, maxfc
                Ref{Cint}, Ref{Cint}, Ref{Cint},                      # trtype, iprint, ncomp
                Ref{Cdouble}, Ptr{Cdouble}, Ref{Cdouble},             # f, g, gpsupn
                Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint},           # iter, fcnt, gcnt, cgcnt
                Ref{Cint}, Ptr{Cint}, Ptr{Cdouble},                   # inform, wi, wd
                Ref{Cdouble},                                         # delmin
            ),
            Cint(n), x, l, u,
            Cint(0), Float64[], Float64[],
            epsgpsn, Cint(nitmax), Cint(nfevalmax),
            Cint(1), Cint(0), Cint(50),
            f, g, gpsupn,
            iter, fcnt, gcnt, cgcnt,
            inform, wi, wd,
            2.0,
        )
        return (; f=f[], gpsupn=gpsupn[], iter=Int(iter[]), fcnt=Int(fcnt[]), gcnt=Int(gcnt[]),
            cgcnt=Int(cgcnt[]), inform=Int(inform[]))
    finally
        CURRENT_FG[] = nothing
        CURRENT_ON_EVAL[] = nothing
    end
end

end # module
