#
# Julia-facing API for the GENCAN FFI wrapper: the @cfunction callback
# bridge, and the public gencan!() entry point (a near drop-in alternative
# to SPGBox.spgbox! at the one call site in packmol_main.jl).
#

# Holds the C function pointer computed in Packmol.__init__() (see
# src/Packmol.jl) — a @cfunction pointer is a JIT address only valid for
# the current process, so it cannot be a top-level `const` value itself.
const GENCAN_FG_CALLBACK = Ref{Ptr{Cvoid}}(C_NULL)

# Fortran callbacks (evalal/evalnal) have no user-data slot, and GENCAN
# calls are always sequential from Julia (never reentrant) — so a single
# global context is safe and is the standard pattern for this kind of
# embedding (same as used by NLopt.jl/Ipopt.jl-style wrappers).
const GENCAN_CURRENT_FG = Ref{Any}(nothing)

# Optional zero-argument hook invoked once per fg! evaluation (i.e. once
# per evalal/evalnal call), for progress reporting. GENCAN has no native
# mid-run callback the way SPGBox does, so this is our own equivalent,
# ticked from inside the FFI bridge callback itself.
const GENCAN_CURRENT_CALLBACK = Ref{Any}(nothing)

#
# The actual C callback body. Must be a plain top-level (non-closure)
# function so the @cfunction pointer computed once in __init__() stays
# valid and reusable across all gencan!() calls.
#
function gencan_fg_callback(n::Cint, x_ptr::Ptr{Cdouble}, f_ptr::Ptr{Cdouble}, g_ptr::Ptr{Cdouble})::Cint
    fg! = GENCAN_CURRENT_FG[]
    try
        nn = Int(n)
        x = unsafe_wrap(Array, x_ptr, nn; own=false)
        g = unsafe_wrap(Array, g_ptr, nn; own=false)
        f = fg!(g, x)
        unsafe_store!(f_ptr, f)
        on_eval = GENCAN_CURRENT_CALLBACK[]
        on_eval === nothing || on_eval()
        return Cint(0)
    catch e
        @error "error in gencan!'s fg! callback" exception = (e, catch_backtrace())
        return Cint(1)
    end
end

struct GencanResult
    f::Float64
    gpsupn::Float64
    iter::Int
    fcnt::Int
    gcnt::Int
    cgcnt::Int
    inform::Int
end

#
# gencan!(fg!, x; nitmax, nfevalmax, epsgpsn, callback) -> GencanResult
#
# Runs the (unconstrained, box bounds ±1e20 — matching how Fortran Packmol
# itself calls GENCAN when no rotation constraints are set) GENCAN
# optimizer on `fg!`, a combined function+gradient closure with the same
# `(g, x) -> f::Float64` signature SPGBox.spgbox! expects. Updates `x`
# in-place with the solution.
#
# `callback`, if given, is a zero-argument function invoked once per fg!
# evaluation (i.e. once per evalal/evalnal call) — used for progress
# reporting. Unlike SPGBox's callback, its return value is ignored: GENCAN
# has no mechanism to stop mid-run from a Julia-side signal, so this can't
# be used for early-exit-on-convergence the way SPGBox's callback is.
#
# GENCAN is hardcoded double precision; `x` must be a Vector{Float64}.
# This is an experimental, local-build-only integration: the first call
# compiles a small shared library from vendored Fortran sources (see
# src/gencan/build.jl), which requires `gfortran` on PATH.
#
function gencan!(
    fg!::F, x::Vector{Float64};
    nitmax::Int=100, nfevalmax::Int=1000, epsgpsn::Float64=1e-6,
    callback::C=nothing,
) where {F,C}
    n = length(x)
    GENCAN_CURRENT_FG[] === nothing || error(
        "gencan!() called while another gencan!() call is already in " *
        "progress — GENCAN's Fortran callbacks have no reentrancy support."
    )
    GENCAN_CURRENT_FG[] = fg!
    GENCAN_CURRENT_CALLBACK[] = callback
    try
        ccall_set_fg_callback(GENCAN_FG_CALLBACK[])

        l = fill(-1.0e20, n)
        u = fill(1.0e20, n)
        g = zeros(n)
        wi = zeros(Int32, n)
        wd = zeros(n * 8)
        delmin = 2.0

        r = ccall_easygencan!(
            n, x, l, u, epsgpsn, nitmax, nfevalmax, g, wi, wd, delmin,
        )
        return GencanResult(r.f, r.gpsupn, r.iter, r.fcnt, r.gcnt, r.cgcnt, r.inform)
    finally
        GENCAN_CURRENT_FG[] = nothing
        GENCAN_CURRENT_CALLBACK[] = nothing
    end
end

@testitem "gencan FFI quadratic" begin
    using Packmol

    if Packmol.gencan_gfortran_available()
        x0 = [1.0, -2.0, 3.0]
        x = zeros(3)
        fg!(g, x) = (g .= 2 .* (x .- x0); sum(abs2, x .- x0))
        result = Packmol.gencan!(fg!, x)
        @test x ≈ x0 atol = 1e-4
        @test result.f < 1e-6
    else
        @warn "gfortran not available: skipping gencan FFI test"
        @test_skip false
    end
end
