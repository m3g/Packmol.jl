#
# Large examples (mostly the tutorials): run on request only, e.g. before a release. See
# README.md. Usage, from the root of the repository:
#
#     julia --project=dev/large_tests -t 14 dev/large_tests/run.jl [options] [case ...]
#
# With no case, all of them are run, one after the other. Options:
#
#     --list          list the cases and their approximate times, and exit
#     --seed N        add `seed N` to the inputs (default: the seed of the input, if any)
#     --workdir DIR   where the inputs are copied and packed (default: a temporary directory)
#
using Packmol
using LazyArtifacts: @artifact_str
using CellListMap: ParticleSystem, pairwise!
using StaticArrays: SVector
using Printf: @sprintf
using Dates: now, format

const root = normpath(joinpath(@__DIR__, "..", ".."))
const tutorials = joinpath(root, "docs", "src", "tutorials", "files")

# name, directory of the input (a function: the artifacts are downloaded only when needed),
# input file (or input files, packed one after the other in the same directory, the last
# one checked and reported, the time being that of all of them), approximate time
# (14 threads)
const cases = [
    ("protein_glycerol", () -> joinpath(tutorials, "protein_glycerol"), "protein_glycerol.inp", "< 1 s"),
    ("water_drop", () -> joinpath(tutorials, "water_drop"), "drop.inp", "< 1 s"),
    ("gaussian_drop", () -> joinpath(tutorials, "water_drop"), "gaussian_drop.inp", "< 1 s"),
    ("zeolite", () -> joinpath(tutorials, "zeolite"), "zeolite.inp", "< 1 s"),
    ("zif8", () -> joinpath(tutorials, "zif8"), "zif8.inp", "< 1 s"),
    ("polymer_dmf", () -> joinpath(tutorials, "polymer_dmf"), "polymer_dmf.inp", "< 1 s"),
    ("nanotube", () -> joinpath(tutorials, "nanotube"), "nanotube.inp", "< 1 s"),
    ("bilayer", () -> joinpath(tutorials, "bilayer"), "bilayer.inp", "2 s"),
    ("nanoliposome_small", () -> joinpath(tutorials, "nanoliposome"), "small.inp", "20 s"),
    ("nanoliposome_large", () -> joinpath(tutorials, "nanoliposome"), "large.inp", "30 s"),
    ("zika", () -> artifact"packmol_zika", "zika.inp", "1 min"),
    ("lnp", () -> artifact"packmol_lnp", "lnp.inp", "5 min"),
    ("exosome", () -> artifact"packmol_exosome", ["proteins.inp", "exosome.inp"], "30 min"),
]

function parse_args(args)
    selected = String[]
    seed = nothing
    workdir = nothing
    list = false
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--list"
            list = true
        elseif a == "--seed"
            seed = parse(Int, args[i+1]); i += 1
        elseif a == "--workdir"
            workdir = args[i+1]; i += 1
        elseif startswith(a, "--")
            error("unknown option $a")
        else
            any(c -> c[1] == a, cases) || error("unknown case $a (see --list)")
            push!(selected, a)
        end
        i += 1
    end
    isempty(selected) && (selected = [c[1] for c in cases])
    workdir = isnothing(workdir) ? mktempdir() : mkpath(workdir)
    return (; selected, seed, workdir, list)
end

# The criterion of the packing, checked independently: no pair of atoms of different
# molecules (not both fixed) closer than the sum of their radii by `tolerance_precision`
# or more. From the positions of the engine (Float64): those of `get_atoms` are Float32,
# with errors of ~1e-4 Å at 1000 Å from the origin, enough to put a pair accepted at
# 1.90002 Å below 1.9. Returns the largest overlap and the minimum distance between
# molecules (molecules are rigid, so the distances within a molecule are those of its
# structure file).
function check_overlaps(sys)
    x = similar(sys.atoms, SVector{3,Float64})
    Packmol.compute_atom_positions!(x, sys.molecule_positions, sys)
    fixed = [sys.structure_types[a.structure_type_index].fixed.fixed for a in sys.atoms]
    cutoff = 2 * maximum(a.radius for a in sys.atoms)
    unitcell = isnothing(sys.unitcell) ? nothing : sys.unitcell
    # (the output, the largest overlap, is reset to zero before the computation; the
    # minimum distance is kept aside, the computation being serial)
    dmin = Ref(cutoff)
    ps = ParticleSystem(; positions = x, cutoff, unitcell, output = 0.0, parallel = false)
    max_overlap = pairwise!(ps) do pair, o
        (; i, j, d) = pair
        ai, aj = sys.atoms[i], sys.atoms[j]
        (ai.molecule_index == aj.molecule_index || (fixed[i] && fixed[j])) && return o
        dmin[] = min(dmin[], d)
        return max(o, ai.radius + aj.radius - d)
    end
    return max_overlap, dmin[]
end

function run_case(name, dir, inputs, workdir, seed)
    dir = dir()
    inputs = inputs isa AbstractString ? [inputs] : inputs
    all(input -> isfile(joinpath(dir, input)), inputs) || return nothing
    work = joinpath(workdir, name)
    rm(work; recursive = true, force = true)
    cp(dir, work)
    chmod(work, 0o755; recursive = true)
    sys, time, log = nothing, 0.0, ""
    for input in inputs
        inp = joinpath(work, input)
        if !isnothing(seed)
            lines = filter(l -> !occursin(r"^\s*seed\b", l), readlines(inp))
            write(inp, "seed $seed\n" * join(lines, '\n') * '\n')
        end
        log = joinpath(work, replace(input, ".inp" => ".log"))
        sys, t = open(log, "w") do io
            redirect_stdout(io) do
                redirect_stderr(io) do
                    t = @elapsed s = Packmol.packmol(inp)
                    (s, t)
                end
            end
        end
        time += t
        sys.status == :packing_ready || break
    end
    s = read(log, String)
    loops = count("Starting packing loop", s)
    evaluations = sum(m -> parse(Int, m[1]), eachmatch(r"Function evaluations:\s+(\d+)", s); init = 0)
    max_overlap, dmin = check_overlaps(sys)
    ok = sys.status == :packing_ready && max_overlap < sys.tolerance_precision
    return (; ok, status = sys.status, loops, evaluations, time, dmin, max_overlap, log)
end

# A small packing, so that the time of the first case does not include the compilation
function warm_up()
    pdb = tempname() * ".pdb"
    write(pdb, "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\nEND\n")
    sys = PackmolSystem([structure_type(pdb; number = 10, constraints = [InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])])];
        output = "", tolerance = 2.0)
    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            Packmol.packmol(sys)
        end
    end
    rm(pdb; force = true)
    return nothing
end

function main(args)
    opts = parse_args(args)
    if opts.list
        for (name, _, input, t) in cases
            println(rpad(name, 22), rpad(input isa AbstractString ? input : join(input, ", "), 30), t)
        end
        return true
    end
    commit = try
        c = readchomp(`git -C $root rev-parse --short HEAD`)
        isempty(readchomp(`git -C $root status --porcelain --untracked-files=no`)) ? c : c * "+dirty"
    catch
        "unknown"
    end
    date = format(now(), "yyyy-mm-dd HH:MM")
    results_file = joinpath(@__DIR__, "results.tsv")
    if !isfile(results_file)
        write(results_file, "date\tcommit\tjulia\tthreads\tcase\tseed\tstatus\tloops\tevaluations\ttime_s\tdmin\n")
    end
    println("Packmol.jl large tests: commit $commit, $(Threads.nthreads()) threads, in $(opts.workdir)")
    warm_up()
    all_ok = true
    for (name, dir, input, t) in cases
        name in opts.selected || continue
        println("  ", rpad(name, 22), "(about $t) ..."); flush(stdout)
        r = run_case(name, dir, input, opts.workdir, opts.seed)
        if isnothing(r)
            println("    skipped: input not found")
            continue
        end
        all_ok &= r.ok
        println(@sprintf("    %s: %s, %d loops, %d evaluations, %.1f s, minimum distance %.5f, largest overlap %.5f",
            r.ok ? "ok" : "FAILED", r.status, r.loops, r.evaluations, r.time, r.dmin, r.max_overlap))
        r.ok || println("    log: ", r.log)
        open(results_file, "a") do io
            println(io, join((date, commit, VERSION, Threads.nthreads(), name, something(opts.seed, ""),
                r.status, r.loops, r.evaluations, round(r.time; digits = 1), round(r.dmin; digits = 3)), '\t'))
        end
    end
    println(all_ok ? "All passed." : "Some cases FAILED.")
    return all_ok
end

main(ARGS) || exit(1)
