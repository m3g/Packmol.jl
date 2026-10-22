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
