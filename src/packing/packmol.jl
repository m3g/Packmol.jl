#
# The `packmol` functions of the API: from an input file, from a file chosen in a file
# dialog, and from a `PackmolSystem` (which runs the packing loop of packing_loop.jl).
#
#
# The keyword arguments of `packmol`, explicit in its signatures (so that they are
# completed in the REPL, and a misspelled one is an error): the global options of the
# input files, set in the system, and the options of the run, passed to `_packmol`. The
# methods for an input file and for a file chosen in the file dialog forward them to
# `packmol(::PackmolSystem)`, and are generated below from these lists.
#
const _packmol_global_options = (
    :(output::Union{Nothing,AbstractString} = nothing),
    :(filetype::Union{Nothing,AbstractString,Symbol} = nothing),
    :(tolerance::Union{Nothing,Real} = nothing),
    :(radscale::Union{Nothing,Real} = nothing),
    :(discale::Union{Nothing,Real} = nothing),
    :(radscale_decay_loops::Union{Nothing,Integer} = nothing),
    :(writeout::Union{Nothing,Real} = nothing),
    :(tolerance_precision::Union{Nothing,Real} = nothing),
    :(constraint_precision::Union{Nothing,Real} = nothing),
    :(maxit::Union{Nothing,Integer} = nothing),
    :(nloop::Union{Nothing,Integer} = nothing),
    :(movefrac::Union{Nothing,Real} = nothing),
    :(nthreads::Union{Nothing,Integer} = nothing),
    :(seed::Union{Nothing,Integer} = nothing),
    :(max_random_init::Union{Nothing,Integer} = nothing),
    :(adjust_constraints_on_init::Union{Nothing,Bool} = nothing),
    :(add_amber_ter::Union{Nothing,Bool} = nothing),
    :(amber_ter_preserve::Union{Nothing,Bool} = nothing),
    :(add_box_sides::Union{Nothing,Bool} = nothing),
    :(periodic_boundary_style::Union{Nothing,AbstractString,Symbol} = nothing),
    :(connect::Union{Nothing,Bool} = nothing),
    :(avoid_overlap::Union{Nothing,Bool} = nothing),
    :(screen_log::Union{Nothing,Bool} = nothing),
    :(check_initial_placement::Union{Nothing,Bool} = nothing),
    :(continue_failed_initial_packing::Union{Nothing,Bool} = nothing),
    :(restart_from::Union{Nothing,AbstractString} = nothing),
    :(pbc = nothing),
    :(unitcell = nothing),
)
const _packmol_run_options = (
    :(active_set::Bool = true),
    :(active_radius::Union{Nothing,Real} = nothing),
    :(active_max_fraction::Real = 0.7),
    :(rotation_scaling::Bool = true),
    :(movebad_tolerance::Real = 0.10),
    :(parallel::Bool = true),
)
# The name of the keyword argument of an expression of the lists above
_kwarg_name(ex::Expr) = ex.args[1] isa Symbol ? ex.args[1] : ex.args[1].args[1]
const _packmol_kwarg_names = map(_kwarg_name, (_packmol_global_options..., _packmol_run_options...))
_kwarg_signature(ex::Expr) = Expr(:kw, ex.args[1], ex.args[2])
# The keyword arguments forwarded by the `packmol` methods of the recipes (see
# ../Recipes/), which set the output, the periodic cell and its shape themselves
const _packmol_recipe_kwargs = filter(ex -> _kwarg_name(ex) ∉ (:output, :pbc, :unitcell),
    (_packmol_global_options..., _packmol_run_options...))
const _packmol_recipe_kwarg_names = map(_kwarg_name, _packmol_recipe_kwargs)

@eval function packmol(input_file::String; D::Int=3, T::DataType=Float64,
    $(map(_kwarg_signature, _packmol_global_options)...),
    $(map(_kwarg_signature, _packmol_run_options)...),
)
    packmol_system = read_packmol_input(input_file; D, T)
    return packmol(packmol_system; $(_packmol_kwarg_names...))
end
@eval function packmol(;
    $(map(_kwarg_signature, _packmol_global_options)...),
    $(map(_kwarg_signature, _packmol_run_options)...),
)
    return packmol(pick_file(); $(_packmol_kwarg_names...))
end

"""
    packmol(input_file::String; options...)

Read a Packmol input file, run the packing optimization, and write the output.
Returns the built `PackmolSystem` (see `packmol(::PackmolSystem)` for details, and for the
keyword arguments, which are the same).

Paths in the input file (structure files, `output`, `restart_from`)
are relative to the directory of the input file.

With no input file, a file explorer is opened to choose it (see [`Packmol.pick_file`](@ref)).
To disable the file explorer, set the environment variable `PACKMOL_GUI="false"`.
"""
packmol(::String)

"""
    packmol(; options...)

Choose a Packmol input file in a file browser, and run it.

The same as `packmol(input_file; options...)`, with the file chosen in the file browser of
[`Packmol.pick_file`](@ref) (disabled by setting the environment variable
`PACKMOL_GUI="false"`, in which case the file is asked for in the terminal).
"""
packmol()

@eval function packmol(packmol_system::PackmolSystem;
    $(map(_kwarg_signature, _packmol_global_options)...),
    $(map(_kwarg_signature, _packmol_run_options)...),
)
    # The global options of the input files, set in the system
    _set_global_options!(packmol_system;
        output, filetype, tolerance, radscale, discale, radscale_decay_loops, writeout,
        tolerance_precision, constraint_precision, maxit, nloop, movefrac, nthreads, seed,
        max_random_init, adjust_constraints_on_init, add_amber_ter, amber_ter_preserve,
        add_box_sides, periodic_boundary_style, connect, avoid_overlap, screen_log, check_initial_placement,
        continue_failed_initial_packing, restart_from, pbc, unitcell,
    )
    ntasks = _ntasks[]
    _ntasks[] = something(packmol_system.nthreads, 0)
    screen_log_before = _screen_log[]
    _screen_log[] = packmol_system.screen_log
    try
        return _packmol(packmol_system;
            active_set, active_max_fraction, rotation_scaling, parallel,
            active_radius = something(active_radius, 2 * packmol_system.tolerance),
            movebad_tolerance,
        )
    finally
        _ntasks[] = ntasks
        _screen_log[] = screen_log_before
    end
end
"""
    packmol(packmol_system::PackmolSystem; options...)

Run the packing optimization on a `PackmolSystem`, mutating it in place (final
`molecule_positions`, ...) and returning that same object — not a bare `Bool` —
so callers can keep using it afterwards (e.g. to postprocess or re-`write_output`
the packed coordinates, or to convert their representation, as with
[`triclinic_to_dodecahedral`](@ref)).

The packing outcome is recorded in `packmol_system.status`: `:packing_ready` if
it converged within `nloop` loops (both `tolerance_precision` and
`constraint_precision` satisfied), `:failed_packing` if it did not, or
`:failed_initial_approximation` if the initial approximation could not satisfy the
constraints (written to the output file; the packing continues only if
`packmol_system.continue_failed_initial_packing` is set), or
`:not_packed` if `packmol_system.check_initial_placement` was set (initial approximation written,
no packing attempted).

**Global options.** Every global keyword of the input files (see [Global keywords](@ref)) can
be given as a keyword argument, with the name of the keyword: `tolerance`, `output`,
`filetype`, `seed`, `nloop`, `maxit`, `movefrac`, `radscale`, `radscale_decay_loops`,
`writeout`, `tolerance_precision`, `constraint_precision`, `nthreads`, `avoid_overlap`, `screen_log`,
`connect`, `add_box_sides`, `add_amber_ter`, `amber_ter_preserve`, `periodic_boundary_style`,
`max_random_init`, `adjust_constraints_on_init`, `check_initial_placement`,
`continue_failed_initial_packing`, `restart_from`, `pbc` and `unitcell`. They are set in
`packmol_system`, replacing the values of the input file (or of the constructor), before
the packing. For example, `packmol(system; tolerance = 2.5, maxit = 50, output = "out.pdb")`.
The values are those of the input keywords: numbers, strings or symbols for the options
(`periodic_boundary_style = :octahedral`), `true`/`false` for the flags; `pbc` takes the
values of the input keyword (`pbc = [40, 40, 40]`, `pbc = (:octahedral, 0, 0, 0, 80)`), and
`unitcell` the six cell parameters or a unit cell matrix. A new `tolerance` sets the radius
of the atoms that have the default one (`tolerance / 2`).

The defaults: `nloop` `500`, `maxit` (optimizer iterations per loop) `30`, `movefrac`
(maximum fraction of the molecules of each type moved per loop) `0.05`, `writeout`
(minimum time, in seconds, between two writes of the output during the packing: the best
solution so far is written after a loop that improves it, if the end of the last write,
or the start of the packing, was at least `writeout` seconds before; the final solution is always
written) `60`, `nthreads` (maximum number of threads used) all the threads of the Julia
session, and `radscale_decay_loops` `20`: the number of loops in which the radii of the
atoms, multiplied at the start by the `radscale` of their structure type (see
[`PackmolSystem`](@ref); by default 5.0 for the free molecules of more than 20 atoms, 1.0
for the others), decay exponentially to their values (to within 0.1% of the initial
excess).

**Options of the run.**

`active_set` (default `true`): after the first loop, optimize only the
molecules with overlaps or constraint violations that may block the convergence, and the
molecules with an atom within `active_radius` (default `2 * tolerance`) of them, with the
other molecules frozen. Used in a loop when the atoms of the active molecules are at most
`active_max_fraction` (default `0.7`) of the atoms of the free molecules; otherwise the
whole system is optimized. Each loop reports the fractions of active molecules and atoms.
Short loops with the active subsystem pack large systems fastest: about 2.5x faster
than whole-system loops of 800 iterations for a nanoliposome of 600,000 atoms.

`rotation_scaling` (default `true`): the optimizer works on rotation variables scaled by
the inertia of each molecule (a unit change moves its atoms by 1 Å RMS), instead of its
Euler angles. With `false`, the Euler angles themselves, much slower for systems of
elongated molecules (lipids, polymers); for comparisons only.

`movebad_tolerance` (default `0.10`): after the first loop, the bad molecules are moved
(see `movefrac`) when a loop improves the objective function by less than
`100 * movebad_tolerance` percent.

`parallel` (default `true`): `false` computes the objective function in a single thread.
"""
packmol(::PackmolSystem)
