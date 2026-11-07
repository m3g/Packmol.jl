#
# Output of the packing loop: the terminal report (progress bar, overlaps) and the
# `<output>_FORCED` file written when the packing fails.
#

#
# Aesthetic line constants (matching Fortran Packmol)
#
const dash_line = repeat('-', 80)
const hash_line = repeat('#', 80)
const _packmol_version = string(pkgversion(@__MODULE__))

#
# When packing fails, write (as Fortran Packmol) the "forced" solution to
# `<output>_FORCED`: the best solution, after an adjustment that satisfies the
# constraints regardless of overlaps. The best solution itself is kept in the
# system (and in the regular output and restart files).
#
function _write_forced_output(packmol_system::PackmolSystem, free_mol_indices::Vector{Int}, RNG, buffers::MemoryBuffers, nloop::Int)
    best_positions = copy(packmol_system.molecule_positions)
    _println()
    _println("  Forcing the solution to fit the constraints...")
    cm_min, cm_max = compute_cm_bounds(packmol_system)
    adjust_constraints!(packmol_system, free_mol_indices, RNG;
        nloop, cm_lo_type=cm_min, cm_hi_type=cm_max, buffers,
    )
    forced_file = packmol_system.output_file * "_FORCED"
    write_output(packmol_system; output_file=forced_file)
    _println("  The forced point was written to the output file: ", forced_file)
    copyto!(packmol_system.molecule_positions, best_positions)
    return forced_file
end

#
# SPGBox callback: progress bar, and stop on convergence
#
function packmol_callback(
    cl_system, tol, iprint, tolerance_precision, constraint_precision, progress_meter;
    f::Union{Nothing,Real} = nothing,
    nfeval::Union{Nothing,Integer} = nothing,
    nfevalmax::Union{Nothing,Integer} = nothing,
    gnorm::Union{Nothing,Real} = nothing,
)
    dmin = min(cl_system.fg.dmin, cl_system.cutoff)
    max_const = cl_system.fg.max_constraint_penalty
    # nfeval/nfevalmax and gnorm are diagnostic only: the progress bar and the
    # values above only update once SPGBox accepts an iteration, so nfeval shows
    # a slow line search, and a blown-up or NaN gnorm a numerical problem.
    next!(progress_meter; showvalues = (
        ("Function value", isnothing(f) ? "" : @sprintf("%.5e", f)),
        ("Minimum distance", _format_dmin(dmin, cl_system.cutoff)),
        ("Maximum constraint violation", @sprintf("%.5e", max_const)),
        ("Function evaluations", isnothing(nfeval) ? "" : "$nfeval / $nfevalmax"),
        ("Projected gradient norm", isnothing(gnorm) ? "" : @sprintf("%.5e", gnorm)),
    ))
    tol_ok = cl_system.fg.max_overlap < tolerance_precision
    const_ok = max_const < constraint_precision
    return tol_ok && const_ok
end

#
# Report of the pairs of atoms overlapping by more than tolerance_precision (the
# convergence criterion), for the pairs of structure types with the most of them
#
function _print_overlaps(packmol_system::PackmolSystem, overlaps::Matrix{Int}; nshow::Int=6)
    total = sum(overlaps)
    total == 0 && return false
    @_printf("  Overlapping pairs of atoms: %d, between molecules of types:\n", total)
    pairs = Tuple{Int,Int,Int}[]
    for j in axes(overlaps, 2), i in 1:j
        overlaps[i, j] > 0 && push!(pairs, (overlaps[i, j], i, j))
    end
    sort!(pairs; rev=true)
    for k in 1:min(nshow, length(pairs))
        n, i, j = pairs[k]
        @_printf("    %10d  %s - %s\n", n, _type_label(packmol_system, i), _type_label(packmol_system, j))
    end
    if length(pairs) > nshow
        @_printf("    %10d  (other pairs of types)\n", sum(first(p) for p in pairs[nshow+1:end]))
    end
    return true
end
function _type_label(packmol_system::PackmolSystem, itype::Int)
    st = packmol_system.structure_types[itype]
    label = string(itype, ": ", basename(st.filename))
    return st.fixed.fixed ? label * " (fixed)" : label
end
