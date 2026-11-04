#
# Packing with other local solvers in place of SPGBox, inside the real packing loop
# (same movebad, stall detection and convergence criteria): the solver of each chunk is
# chosen by redefining `Packmol._local_solver`.
#
#   julia --project=dev/solvers -t 14 dev/solvers/run.jl <solver> <input.inp> [seed]
#
# <solver>: spgbox, gencan, lbfgs. Environment variables: ROTATION_SCALING (default true),
# LBFGS_M (memory, default 10), LBFGS_LS, LBFGS_ALPHA, GENCAN_MAXIT (iterations per chunk,
# default 20, as in Fortran Packmol), NLOOP, and GRADCHECK or LSTRACE (diagnostics of the
# first chunk).
#
# Prints, at the end, the total number of objective evaluations of the packing loop
# (function and gradient are always evaluated together), the time and the status.
#
using Packmol, SPGBox, Optim, NLSolversBase, LineSearches
include(joinpath(@__DIR__, "gencan", "GencanFFI.jl"))
using .GencanFFI

const SOLVER = Symbol(ARGS[1])
const INPUT = ARGS[2]
const SEED = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : nothing
const LBFGS_M = parse(Int, get(ENV, "LBFGS_M", "10"))
const GENCAN_MAXIT = parse(Int, get(ENV, "GENCAN_MAXIT", "20"))
# L-BFGS line search: hz (Hager-Zhang, Optim's default), bt (backtracking, cubic), mt (More-Thuente)
const LBFGS_LS = get(ENV, "LBFGS_LS", "hz")
# Length of the first trial step of each line search (scaled: the step of the first
# iteration has this norm), or the unit step of Optim's default when not set
const LBFGS_ALPHA = haskey(ENV, "LBFGS_ALPHA") ? parse(Float64, ENV["LBFGS_ALPHA"]) : nothing
# Scaling of the rotations of the packing loop (src/rigid_body/rotation_metric.jl):
# ROTATION_SCALING=false packs with the Euler angles themselves
const ROTATION_SCALING = get(ENV, "ROTATION_SCALING", "true") == "true"
# Number of packing loops (default: the input's)
const NLOOP = haskey(ENV, "NLOOP") ? parse(Int, ENV["NLOOP"]) : nothing
# Iterations of the solver, over all chunks
const NIT = Ref(0)
SOLVER in (:spgbox, :gencan, :lbfgs) || error("unknown solver $SOLVER")

# Objective evaluations of the packing loop, over all chunks
const NEVAL = Ref(0)
# Chunks so far (GRADCHECK=n checks the gradient at the start of chunk n)
const NCHUNK = Ref(0)

# Directional derivatives along -g (the steepest descent direction), along random
# directions, and along the center-of-mass and the angle coordinates separately,
# compared with finite differences of f, at the point of the first chunk; then stops.
function gradcheck(fg!, x)
    n = length(x)
    g = similar(x)
    f0 = fg!(g, x)
    g0 = copy(g)
    gt = similar(x)
    cm_mask = [mod1(i, 6) <= 3 for i in 1:n]
    dirs = [
        ("-g", -g0),
        ("-g (cm only)", -g0 .* cm_mask),
        ("-g (angles only)", -g0 .* .!cm_mask),
        ("random", randn(n)),
    ]
    println("GRADCHECK f0 = $f0, |g| = $(sqrt(sum(abs2, g0)))")
    for (name, d) in dirs
        d ./= sqrt(sum(abs2, d))
        slope = sum(g0 .* d)
        println("GRADCHECK direction $name: g.d = $slope")
        for h in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
            fp = fg!(gt, x .+ h .* d)
            fm = fg!(gt, x .- h .* d)
            println("GRADCHECK   h = $h: central = $((fp - fm) / (2h))  forward = $((fp - f0) / h)")
        end
    end
    fg!(g, x)
    exit(0)
end

# Evaluations and step lengths of each L-BFGS iteration from the point of the first
# chunk (60 iterations), with the line searches of LBFGS_LS; then stops.
function lstrace(fg!, x)
    nf = Ref(0)
    gbuf = similar(x)
    fg3! = (F, G, x) -> begin
        nf[] += 1
        f = fg!(G === nothing ? gbuf : G, x)
        F === nothing ? nothing : f
    end
    linesearch = LBFGS_LS == "bt" ? LineSearches.BackTracking(order=3) :
        LBFGS_LS == "mt" ? LineSearches.MoreThuente() : LineSearches.HagerZhang()
    alphaguess = isnothing(LBFGS_ALPHA) ? LineSearches.InitialStatic() :
        LineSearches.InitialStatic(alpha=LBFGS_ALPHA, scaled=true)
    last = Ref(0)
    cb = state -> begin
        println("LSTRACE f=$(state.f_x) |g|=$(sqrt(sum(abs2, state.g_x))) ",
            "alpha=$(state.alpha) |dx|=$(sqrt(sum(abs2, state.dx))) evals=$(nf[] - last[])")
        last[] = nf[]
        false
    end
    options = Optim.Options(; iterations=60, f_calls_limit=100000, g_abstol=1e-8, callback=cb)
    r = optimize(NLSolversBase.only_fg!(fg3!), copy(x), LBFGS(; m=LBFGS_M, linesearch, alphaguess), options)
    println("LSTRACE total evals $(nf[]) iterations $(Optim.iterations(r)) f $(Optim.minimum(r))")
    exit(0)
end

struct ChunkResult
    f::Float64
    ierr::Int     # 0: converged (or stopped by the callback), 1: iterations, 2: evaluations
    nfeval::Int
    gnorm::Float64
end

function Packmol._local_solver(
    fg!, x::Vector{Float64}, lower, upper;
    callback, vaux=nothing, nitmax::Int, nfevalmax::Int,
)
    nf = Ref(0)
    counted = (g, x) -> (NEVAL[] += 1; nf[] += 1; fg!(g, x))
    NCHUNK[] += 1
    haskey(ENV, "GRADCHECK") && NCHUNK[] == parse(Int, ENV["GRADCHECK"]) && gradcheck(fg!, x)
    haskey(ENV, "LSTRACE") && lstrace(fg!, x)
    if SOLVER === :spgbox
        r = SPGBox.spgbox!(counted, x; callback, vaux, nitmax, nfevalmax, lower, upper)
        NIT[] += r.nit
        return r
    elseif SOLVER === :gencan
        xstop = similar(x)
        stopped = Ref(false)
        on_eval = (xe, g, f) -> begin
            stop = callback((; f, nfeval=nf[], gnorm=maximum(abs, g)))::Bool
            stop && (xstop .= xe; stopped[] = true)
            stop
        end
        r = GencanFFI.gencan!(counted, x; lower, upper, nitmax=GENCAN_MAXIT, nfevalmax, on_eval)
        stopped[] && (x .= xstop)
        NIT[] += r.iter
        ierr = (stopped[] || r.inform in (0, 1)) ? 0 : r.inform == 8 ? 2 : 1
        return ChunkResult(r.f, ierr, nf[], r.gpsupn)
    else # lbfgs: unbounded (the centers of mass are wrapped into the cell after every loop)
        gbuf = similar(x)
        fg3! = (F, G, x) -> begin
            f = counted(G === nothing ? gbuf : G, x)
            F === nothing ? nothing : f
        end
        stopped = Ref(false)
        cb = state -> begin
            stop = callback((; f=state.f_x, nfeval=nf[], gnorm=maximum(abs, state.g_x)))::Bool
            stopped[] |= stop
            stop
        end
        options = Optim.Options(; iterations=nitmax, f_calls_limit=nfevalmax, g_abstol=1e-6, callback=cb)
        linesearch = LBFGS_LS == "bt" ? LineSearches.BackTracking(order=3) :
            LBFGS_LS == "mt" ? LineSearches.MoreThuente() : LineSearches.HagerZhang()
        alphaguess = isnothing(LBFGS_ALPHA) ? LineSearches.InitialStatic() :
            LineSearches.InitialStatic(alpha=LBFGS_ALPHA, scaled=true)
        r = optimize(NLSolversBase.only_fg!(fg3!), x, LBFGS(; m=LBFGS_M, linesearch, alphaguess), options)
        NIT[] += Optim.iterations(r)
        x .= Optim.minimizer(r)
        ierr = (stopped[] || Optim.converged(r)) ? 0 : nf[] >= nfevalmax ? 2 : 1
        return ChunkResult(Optim.minimum(r), ierr, nf[], 0.0)
    end
end

# The seed is appended to a copy of the input, next to it (its paths are relative to it)
input = INPUT
if !isnothing(SEED)
    input = joinpath(dirname(abspath(INPUT)), "_solver_seed$(SEED).inp")
    write(input, read(INPUT, String) * "\nseed $SEED\n" * (isnothing(NLOOP) ? "" : "nloop $NLOOP\n"))
end
t0 = time()
system = packmol(input; rotation_scaling=ROTATION_SCALING)
println("RESULT solver=$SOLVER rotation_scaling=$ROTATION_SCALING m=$LBFGS_M ls=$LBFGS_LS alpha=$LBFGS_ALPHA gencan_maxit=$GENCAN_MAXIT seed=$SEED nit=$(NIT[]) ",
    "time=$(round(time() - t0; digits=1)) s evals=$(NEVAL[]) status=$(system.status)")
