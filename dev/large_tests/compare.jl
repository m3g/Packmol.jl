#
# Compares the results of two versions of Packmol.jl in the large tests: the runs of
# run.jl with several seeds, one set of runs per version, labeled. For a change with
# widespread consequences (e.g. the heuristics of the packing), from the root of the
# repository:
#
#     git checkout main
#     julia --project=dev/large_tests -t 14 dev/large_tests/run.jl --seeds 1:3 --skip exosome --label main --results validation.tsv
#     git checkout <branch>
#     julia --project=dev/large_tests -t 14 dev/large_tests/run.jl --seeds 1:3 --skip exosome --label new --results validation.tsv
#     julia dev/large_tests/compare.jl validation.tsv main new
#
# For each case, the runs that converged (and passed the independent check of run.jl),
# and the medians of the loops, function evaluations and times of the runs. A case is
# flagged when it converged in fewer runs with the second version, or when its median
# time grew by more than 30% (and 2 s).
#
using Printf: @sprintf

function median(x)
    v = sort(collect(x))
    n = length(v)
    return isodd(n) ? float(v[(n+1)÷2]) : (v[n÷2] + v[n÷2+1]) / 2
end

function read_results(file)
    lines = readlines(file)
    header = split(lines[1], '\t')
    col(name) = findfirst(==(name), header)
    rows = [split(l, '\t') for l in lines[2:end] if !isempty(strip(l))]
    return [(; label = r[col("label")], case = r[col("case")], status = r[col("status")],
        loops = parse(Int, r[col("loops")]), evaluations = parse(Int, r[col("evaluations")]),
        time = parse(Float64, r[col("time_s")])) for r in rows]
end

function summary(rows, label, case)
    rs = filter(r -> r.label == label && r.case == case, rows)
    isempty(rs) && return nothing
    ok = count(r -> r.status == "packing_ready", rs)
    return (; n = length(rs), ok, loops = median(r.loops for r in rs),
        evaluations = median(r.evaluations for r in rs), time = median(r.time for r in rs))
end

function main(args)
    length(args) == 3 || error("usage: julia compare.jl <results.tsv> <label 1> <label 2>")
    file, a, b = args
    rows = read_results(file)
    cases = unique(r.case for r in rows if r.label in (a, b))
    fmt(s) = isnothing(s) ? rpad("-", 34) :
        rpad(@sprintf("%d/%d %6.0f %9.0f %8.1f", s.ok, s.n, s.loops, s.evaluations, s.time), 34)
    println(rpad("case", 24), rpad("$a: ok loops evals time(s)", 34), rpad("$b: ok loops evals time(s)", 34))
    nflagged = 0
    for case in cases
        sa = summary(rows, a, case)
        sb = summary(rows, b, case)
        flag = ""
        if !isnothing(sa) && !isnothing(sb)
            if sb.ok / sb.n < sa.ok / sa.n
                flag = "  <- converged less"
            elseif sb.time > 1.3 * sa.time && sb.time - sa.time > 2
                flag = "  <- slower"
            end
        end
        isempty(flag) || (nflagged += 1)
        println(rpad(case, 24), fmt(sa), fmt(sb), flag)
    end
    println(nflagged == 0 ? "No regressions." : "$nflagged case(s) flagged.")
    return nflagged == 0
end

main(ARGS) || exit(1)
