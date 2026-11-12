#
# Output of the packing loop: the terminal report (progress bar, overlaps) and the
# `<name>-FORCED.<ext>` file written when the packing fails.
#

#
# Aesthetic line constants (matching Fortran Packmol)
#
const dash_line = repeat('-', 80)
const hash_line = repeat('#', 80)
const _packmol_version = string(pkgversion(@__MODULE__))

#
# When packing fails, write (as Fortran Packmol) the "forced" solution to
# `<name>-FORCED.<ext>` (`_forced_file`; Fortran: `<output>_FORCED`): the best solution, after an adjustment that satisfies the
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
    forced_file = _forced_file(packmol_system.output_file)
    write_output(packmol_system; output_file=forced_file)
    _println("  The forced point was written to the output file: ", forced_file)
    copyto!(packmol_system.molecule_positions, best_positions)
    return forced_file
end

# `system.pdb` -> `system-FORCED.pdb`: the extension kept, so the format is that of the output
function _forced_file(output_file::AbstractString)
    base, ext = splitext(output_file)
    return string(base, "-FORCED", ext)
end

#
# SPGBox callback: progress bar, and stop on convergence
#
function packmol_callback(cl_system, tolerance_precision, constraint_precision, progress_meter)
    next!(progress_meter)
    tol_ok = cl_system.fg.max_overlap < tolerance_precision
    const_ok = cl_system.fg.max_constraint_penalty < constraint_precision
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
