#
# `@_sync`: as `Base.@sync`, for the tasks created with `Threads.@spawn`, but
# collecting them in a `Vector{Task}` instead of a `Channel{Any}`, so that
# waiting for them compiles with `juliac --trim` (`Base.@sync` waits on the
# untyped elements of its channel). `@spawn` registers each task in the
# enclosing `@sync` by `put!`-ing it into the local variable named
# `Base.sync_varname`, which is defined here as a `_SyncTasks` instead.
#
struct _SyncTasks
    tasks::Vector{Task}
end
Base.put!(s::_SyncTasks, task::Task) = push!(s.tasks, task)

macro _sync(block)
    var = esc(Base.sync_varname)
    quote
        let $var = _SyncTasks(Task[])
            v = $(esc(block))
            for task in $var.tasks
                wait(task)
            end
            v
        end
    end
end

#
# Number of tasks the parallel loops are split into, and hence the maximum number of
# threads that work (the `nthreads` input keyword and keyword argument of `packmol`,
# or the `-t` flag of the executable, which is built with all threads available).
# 0 means all threads.
#
const _ntasks = Ref(0)
_nbatches() = _ntasks[] > 0 ? min(_ntasks[], Threads.nthreads()) : Threads.nthreads()

# `nbatches` of CellListMap's `ParticleSystem`: its heuristic (0, 0) unless the number
# of tasks is limited.
_cl_nbatches() = _nbatches() < Threads.nthreads() ? (_nbatches(), _nbatches()) : (0, 0)
