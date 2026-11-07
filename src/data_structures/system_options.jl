#
# The options of a `PackmolSystem` given as keyword arguments, or by the global keywords of
# an input file: the periodic cell, the tolerance, and the other global options.
#
# The periodic cell of the `pbc` keyword (`values`: its values, as strings)
#   3 values (side lengths): the box from the origin to (a, b, c)
#   6 values (xmin ymin zmin xmax ymax zmax): center at the midpoint
#   dodecahedral|octahedral cx cy cz d
#
function _set_pbc!(sys::PackmolSystem{D,T}, values::AbstractVector{<:AbstractString}) where {D,T}
    if !isempty(values) && values[1] == "dodecahedral"
        sys.unitcell, sys.unitcell_center = parse_pbc_dodecahedral(T, values[2:end], Val(D))
        return sys
    end
    if !isempty(values) && values[1] == "octahedral"
        sys.unitcell, sys.unitcell_center = parse_pbc_octahedral(T, values[2:end], Val(D))
        return sys
    end
    vals = [_parse_value(T, "pbc", v) for v in values]
    if length(vals) == D
        # 3 values: side lengths of the box from the origin to (a, b, c), as in
        # Fortran Packmol (equivalent to `pbc 0 0 0 a b c`)
        sys.unitcell = Matrix{T}(Diagonal(vals))
        sys.unitcell_center = SVector{D,T}(vals) / 2
    elseif length(vals) == 2 * D
        lo = SVector{D,T}(@view(vals[1:D]))
        hi = SVector{D,T}(@view(vals[D+1:2*D]))
        sys.unitcell = Matrix{T}(Diagonal(Vector(hi - lo)))
        sys.unitcell_center = (lo + hi) / 2
    else
        throw(ArgumentError("pbc keyword requires $D (side lengths) or $(2*D) (min/max coords) values, got $(length(vals))"))
    end
    return sys
end

# The periodic cell of the `unitcell` keyword: a b c α β γ (CRYST1 convention), centered
# at the origin
function _set_unitcell!(sys::PackmolSystem{D,T}, values::AbstractVector{<:AbstractString}) where {D,T}
    if D != 3 || length(values) != 6
        throw(ArgumentError("unitcell keyword requires 6 values (a b c α β γ), got $(length(values))"))
    end
    p = [_parse_value(T, "unitcell", v) for v in values]
    sys.unitcell = Matrix{T}(unitcell_matrix(T, p[1], p[2], p[3], p[4], p[5], p[6]))
    sys.unitcell_center = zero(SVector{D,T})
    return sys
end

#
# A new tolerance: the atoms with the default radius (tolerance/2) get the new default
#
function _set_tolerance!(sys::PackmolSystem{D,T}, tolerance::Real) where {D,T}
    old, new = sys.tolerance / 2, T(tolerance) / 2
    for st in sys.structure_types
        replace!(r -> r == old ? new : r, st.radii)
    end
    sys.atoms = [a.radius == old ? AtomData{T}(a.molecule_index, a.structure_type_index, new, a.constraints) : a
                 for a in sys.atoms]
    sys.tolerance = T(tolerance)
    return sys
end

# The values of a `pbc` or `unitcell` keyword argument as those of the input keyword:
# numbers, and the cell shape (`:octahedral`, `"dodecahedral"`) as the first one
_keyword_values(v) = String[string(x) for x in v]

#
# The global options of the input files given as keyword arguments of
# `packmol(::PackmolSystem)`: they are set in the system (see `packmol`).
#
function _set_global_options!(sys::PackmolSystem{D,T};
    output::Union{Nothing,AbstractString} = nothing,
    filetype::Union{Nothing,AbstractString,Symbol} = nothing,
    tolerance::Union{Nothing,Real} = nothing,
    radscale::Union{Nothing,Real} = nothing,
    discale::Union{Nothing,Real} = nothing,
    radscale_decay_loops::Union{Nothing,Integer} = nothing,
    writeout::Union{Nothing,Real} = nothing,
    tolerance_precision::Union{Nothing,Real} = nothing,
    constraint_precision::Union{Nothing,Real} = nothing,
    maxit::Union{Nothing,Integer} = nothing,
    nloop::Union{Nothing,Integer} = nothing,
    movefrac::Union{Nothing,Real} = nothing,
    nthreads::Union{Nothing,Integer} = nothing,
    seed::Union{Nothing,Integer} = nothing,
    max_random_init::Union{Nothing,Integer} = nothing,
    adjust_constraints_on_init::Union{Nothing,Bool} = nothing,
    add_amber_ter::Union{Nothing,Bool} = nothing,
    amber_ter_preserve::Union{Nothing,Bool} = nothing,
    add_box_sides::Union{Nothing,Bool} = nothing,
    periodic_boundary_style::Union{Nothing,AbstractString,Symbol} = nothing,
    connect::Union{Nothing,Bool} = nothing,
    avoid_overlap::Union{Nothing,Bool} = nothing,
    screen_log::Union{Nothing,Bool} = nothing,
    check_initial_placement::Union{Nothing,Bool} = nothing,
    continue_failed_initial_packing::Union{Nothing,Bool} = nothing,
    restart_from::Union{Nothing,AbstractString} = nothing,
    pbc = nothing,
    unitcell = nothing,
) where {D,T}
    isnothing(output) || (sys.output_file = expanduser(String(output)))
    isnothing(filetype) || _set_input_keyword!(sys, "filetype", string(filetype))
    isnothing(tolerance) || _set_tolerance!(sys, tolerance)
    isnothing(discale) || _set_input_keyword!(sys, "discale", string(discale))
    isnothing(radscale) || (sys.radscale = T(radscale))
    isnothing(radscale_decay_loops) || (sys.radscale_decay_loops = Int(radscale_decay_loops))
    isnothing(writeout) || (sys.writeout = T(writeout))
    isnothing(tolerance_precision) || (sys.tolerance_precision = T(tolerance_precision))
    isnothing(constraint_precision) || (sys.constraint_precision = T(constraint_precision))
    isnothing(maxit) || (sys.max_iter = Int(maxit))
    isnothing(nloop) || (sys.nloop = Int(nloop))
    isnothing(movefrac) || (sys.movefrac = T(movefrac))
    isnothing(nthreads) || _set_input_keyword!(sys, "nthreads", string(nthreads))
    isnothing(seed) || (sys.seed = Int(seed))
    isnothing(max_random_init) || (sys.max_random_init = Int(max_random_init))
    isnothing(adjust_constraints_on_init) || (sys.adjust_constraints_on_init = adjust_constraints_on_init)
    isnothing(add_amber_ter) || (sys.add_amber_ter = add_amber_ter)
    isnothing(amber_ter_preserve) || (sys.amber_ter_preserve = amber_ter_preserve)
    isnothing(add_box_sides) || (sys.add_box_sides = add_box_sides)
    isnothing(periodic_boundary_style) || _set_input_keyword!(sys, "periodic_boundary_style", string(periodic_boundary_style))
    isnothing(connect) || (sys.connect = connect)
    isnothing(avoid_overlap) || (sys.avoid_overlap = avoid_overlap)
    isnothing(screen_log) || (sys.screen_log = screen_log)
    isnothing(check_initial_placement) || (sys.check_initial_placement = check_initial_placement)
    isnothing(continue_failed_initial_packing) || (sys.continue_failed_initial_packing = continue_failed_initial_packing)
    isnothing(restart_from) || (sys.restart_from = expanduser(String(restart_from)))
    if !isnothing(pbc)
        _set_pbc!(sys, _keyword_values(pbc))
    elseif !isnothing(unitcell)
        if unitcell isa AbstractMatrix
            sys.unitcell = Matrix{T}(unitcell)
            sys.unitcell_center = zero(SVector{D,T})
        else
            _set_unitcell!(sys, _keyword_values(unitcell))
        end
    end
    if !isnothing(pbc) || !isnothing(unitcell)
        sys.wigner_seitz_cell = _is_wigner_seitz_cell(something(sys.unitcell))
        _check_fixed_inside_cell(sys)
    end
    isnothing(restart_from) || _check_restart_files(sys)
    return sys
end
