#
# Raw ccall wrapper around the compiled GENCAN shared library.
#
# gfortran's default `integer` kind is 4 bytes for this vendored source
# (no -fdefault-integer-8), so every GENCAN integer argument MUST be
# marshaled as Cint (Int32), not Int/Int64 — using the wrong width here is
# a silent stack-corruption bug, not a type error caught at compile time.
#

function ccall_set_fg_callback(cptr::Ptr{Cvoid})
    fptr = gencan_symbol(:set_fg_callback)
    ccall(fptr, Cvoid, (Ptr{Cvoid},), cptr)
    return nothing
end

#
# Thin wrapper around the Fortran `easygencan` subroutine (gfortran-mangled
# symbol `easygencan_`). All arrays are passed by pointer (matching F77's
# pass-by-reference convention, the same pattern Julia uses to call BLAS/
# LAPACK directly); all scalars are passed as Ref{Cint}/Ref{Cdouble}.
#
# On entry: x is the initial point (n doubles), l/u are box bounds (n
# doubles each, use ±1e20 for "unconstrained"), wi/wd are work arrays
# (n Int32, 8n Float64).
#
# On return: x holds the solution, f the objective value, g the gradient,
# gpsupn the sup-norm of the projected gradient, iter/fcnt/gcnt/cgcnt the
# iteration/eval counters, inform the termination code (see gencan.f
# easygencan docstring: 0/1 = converged, other positive = other stopping
# criterion, negative = error in the user callback).
#
function ccall_easygencan!(
    n::Int, x::Vector{Float64}, l::Vector{Float64}, u::Vector{Float64},
    epsgpsn::Float64, maxit::Int, maxfc::Int,
    g::Vector{Float64}, wi::Vector{Int32}, wd::Vector{Float64}, delmin::Float64,
)
    fptr = gencan_symbol(:easygencan_)

    n_c = Cint(n)
    m_c = Cint(0)
    lambda = Float64[]
    rho = Float64[]
    epsgpsn_c = epsgpsn
    maxit_c = Cint(maxit)
    maxfc_c = Cint(maxfc)
    trtype_c = Cint(1)
    iprint_c = Cint(0)
    ncomp_c = Cint(50)
    f = Ref{Float64}(0.0)
    gpsupn = Ref{Float64}(0.0)
    iter = Ref{Cint}(0)
    fcnt = Ref{Cint}(0)
    gcnt = Ref{Cint}(0)
    cgcnt = Ref{Cint}(0)
    inform = Ref{Cint}(0)
    delmin_c = delmin

    ccall(
        fptr, Cvoid,
        (
            Ref{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},   # n,x,l,u
            Ref{Cint}, Ptr{Cdouble}, Ptr{Cdouble},                  # m,lambda,rho
            Ref{Cdouble}, Ref{Cint}, Ref{Cint},                     # epsgpsn,maxit,maxfc
            Ref{Cint}, Ref{Cint}, Ref{Cint},                        # trtype,iprint,ncomp
            Ref{Cdouble}, Ptr{Cdouble}, Ref{Cdouble},               # f,g,gpsupn
            Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint},             # iter,fcnt,gcnt,cgcnt
            Ref{Cint}, Ptr{Cint}, Ptr{Cdouble},                     # inform,wi,wd
            Ref{Cdouble},                                            # delmin
        ),
        n_c, x, l, u,
        m_c, lambda, rho,
        epsgpsn_c, maxit_c, maxfc_c,
        trtype_c, iprint_c, ncomp_c,
        f, g, gpsupn,
        iter, fcnt, gcnt, cgcnt,
        inform, wi, wd,
        delmin_c,
    )

    return (
        f=f[], gpsupn=gpsupn[], iter=Int(iter[]), fcnt=Int(fcnt[]),
        gcnt=Int(gcnt[]), cgcnt=Int(cgcnt[]), inform=Int(inform[]),
    )
end
