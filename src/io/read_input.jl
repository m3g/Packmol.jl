#
# The input file parser: the global keywords, and the structure blocks (whose contents are
# read by read_structure_block.jl).
#
#=
    _parse_value

Parses the input value for a keyword. If the value cannot be parsed, an error is thrown.

The `_val_check` function is used to check the value after parsing, with different 
optional functions for different values. The function must return an error if the
value is not valid for the given input keyword.

=#
function _parse_value(::Type{S}, keyword::String, input_value::AbstractString; _val_check=x -> nothing) where {S}
    value = tryparse(S, input_value)
    if isnothing(value)
        throw(ArgumentError("""\n
            Could not parse value $input_value for keyword $keyword. Expected type: $S.

        """))
    end
    _val_check(value)
    return value
end
_parse_value(::Type{String}, keyword::String, input_value::AbstractString) = String(input_value)

function _parse_options(name::String, val::AbstractString, options::Tuple{Pair{String,V},Vararg{Pair{String,V}}}) where {V}
    for pair in options
        if first(pair) == val
            return last(pair)
        end
    end
    throw(ArgumentError("""\n
        Value "$val" is not valid for option: $name.

    """))
end

const _yes_no = ("yes" => true, "no" => false)

# Warnings of the input file parser: printed directly, rather than through the
# logging macros, to keep the parser compilable with `juliac --trim`.
_warn(msg::AbstractString) = _screen_log[] ? _with_stderr(io -> println(io, "Warning: ", msg)) : nothing

#
# Set the field of `sys` corresponding to a global keyword of the input file,
# with value `val` (an empty string for the value-less flag keywords, e.g.
# `add_box_sides`). Returns `false` if `keyword` is not one of them. An
# explicit chain of comparisons, rather than a dictionary of closures, so
# that each field is set with a value of its own (inferred) type.
#
function _set_input_keyword!(sys::PackmolSystem{D,T}, keyword::AbstractString, val::AbstractString) where {D,T}
    if keyword == "filetype"
        sys.filetype = _parse_options("filetype", val, ("pdb" => "pdb", "xyz" => "xyz", "mmcif" => "mmcif", "cif" => "mmcif"))
    elseif keyword == "output"
        sys.output_file = expanduser(_parse_value(String, "output", val))
    elseif keyword == "tolerance"
        sys.tolerance = _parse_value(T, "tolerance", val)
    elseif keyword == "radscale"
        sys.radscale = _parse_value(T, "radscale", val)
    elseif keyword == "writeout"
        sys.writeout = _parse_value(T, "writeout", val)
    elseif keyword == "radscale_decay_loops"
        sys.radscale_decay_loops = _parse_value(Int, "radscale_decay_loops", val)
    elseif keyword == "discale"
        _warn("discale is a legacy keyword name; use radscale instead.")
        sys.radscale = _parse_value(T, "discale", val)
    elseif keyword == "tolerance_precision"
        sys.tolerance_precision = _parse_value(T, "tolerance_precision", val)
    elseif keyword == "constraint_precision"
        sys.constraint_precision = _parse_value(T, "constraint_precision", val)
    elseif keyword == "maxit"
        sys.max_iter = _parse_value(Int, "maxit", val)
    elseif keyword == "nloop"
        sys.nloop = _parse_value(Int, "nloop", val)
    elseif keyword == "movefrac"
        sys.movefrac = _parse_value(T, "movefrac", val)
    elseif keyword == "nthreads"
        sys.nthreads = _parse_value(Int, "nthreads", val)
        sys.nthreads >= 1 || throw(ArgumentError("""\n
            nthreads must be at least 1, got $(sys.nthreads).

        """))
    elseif keyword == "max_random_init"
        sys.max_random_init = _parse_value(Int, "max_random_init", val)
    elseif keyword == "adjust_constraints_on_init"
        sys.adjust_constraints_on_init = _parse_options("adjust_constraints_on_init", val, _yes_no)
    elseif keyword == "add_amber_ter"
        sys.add_amber_ter = true
    elseif keyword == "amber_ter_preserve"
        sys.amber_ter_preserve = true
    elseif keyword == "periodic_boundary_style"
        sys.periodic_boundary_style = _parse_options("periodic_boundary_style", val,
            ("triclinic" => :triclinic, "dodecahedral" => :dodecahedral, "octahedral" => :octahedral))
    elseif keyword == "add_box_sides"
        sys.add_box_sides = true
    elseif keyword == "connect"
        sys.connect = _parse_options("connect", val, _yes_no)
    elseif keyword == "ignore_conect"
        sys.connect = false
    elseif keyword == "seed"
        sys.seed = _parse_value(Int, "seed", val)
    elseif keyword == "avoid_overlap"
        sys.avoid_overlap = _parse_options("avoid_overlap", val, _yes_no)
    elseif keyword == "screen_log"
        sys.screen_log = _parse_options("screen_log", val, _yes_no)
    elseif keyword == "check_initial_placement"
        sys.check_initial_placement = true
    elseif keyword == "check"
        _warn("check is the legacy keyword name; use check_initial_placement instead.")
        sys.check_initial_placement = true
    elseif keyword == "continue_failed_initial_packing"
        sys.continue_failed_initial_packing = true
    elseif keyword == "voxel_place"
        sys.voxel_place = isempty(strip(val)) ? 1000 : _parse_value(Int, "voxel_place", val)
        sys.voxel_place >= 1 || throw(ArgumentError("""\n
            voxel_place: max_tries_mult must be at least 1, got $(sys.voxel_place).

        """))
    elseif keyword == "restart_from"
        sys.restart_from = expanduser(_parse_value(String, "restart_from", val))
    else
        return false
    end
    return true
end


# Fortran Packmol keywords that have no use in Packmol.jl: they are accepted,
# so that old input files still run, with a warning that they are ignored.
const packmol_legacy_keywords = Dict{String,String}(
    "fscale" => "fscale is a deprecated keyword and was ignored.",
    "fbins" => "fbins is a deprecated keyword and was ignored.",
    "iprint1" => "iprint1 is a deprecated keyword and was ignored.",
    "iprint2" => "iprint2 is a deprecated keyword and was ignored.",
    "iprint3" => "iprint3 is a deprecated keyword and was ignored.",
    "precision" => "precision is a deprecated keyword and was ignored. Use instead: tolerance_precision and/or constraint_precision",
    "packall" => "packall is a deprecated keyword and was ignored (it is the only option).",
    "use_short_tol" => "use_short_tol is a deprecated keyword and was ignored.",
    "short_tol_dist" => "short_tol_dist is a deprecated keyword and was ignored.",
    "short_tol_scale" => "short_tol_scale is a deprecated keyword and was ignored.",
    "short_radius" => "short_radius is a deprecated keyword and was ignored.",
    "short_radius_scale" => "short_radius_scale is a deprecated keyword and was ignored.",
    "sidemax" => "sidemax is a deprecated keyword and was ignored.",
    "maxmove" => "maxmove is a deprecated keyword and was ignored.",
    "nloop0" => "nloop0 is a deprecated keyword and was ignored.",
    "disable_movebad" => "disable_movebad is a deprecated keyword and was ignored.",
    "writecrd" => "writecrd is a deprecated keyword and was ignored (CHARMM CRD output is not supported).",
    "non_standard_conect" => "non_standard_conect is a deprecated keyword and was ignored (CONECT records are read in the standard PDB format).",
    "hexadecimal_indices" => "hexadecimal_indices is a deprecated keyword and was ignored (atom indices above 99999 are always written in hexadecimal).",
    "movebadrandom" => "movebadrandom is a deprecated keyword and was ignored (bad molecules are always moved to random positions).",
    "restart_to" => "restart_to is a deprecated keyword and was ignored: to restart, give the output file (PDB, XYZ or mmCIF) to restart_from.",
    "randominitialpoint" => "randominitialpoint is a deprecated keyword and was ignored.",
    "optimization_print_level" => "optimization_print_level is a deprecated keyword and was ignored.",
    "chkgrad" => "chkgrad is a deprecated keyword and was ignored.",
    "writebad" => "writebad is a deprecated keyword and was ignored (the output is written after every loop that improves the solution).",
)

# Keywords recognized by Fortran Packmol but not yet wired into Packmol.jl
# (see `dev/drop-in-replacement.md`). Listed here, rather than in
# `_set_input_keyword!`, so that using them produces a clear warning instead
# of a crash while they remain unimplemented.
const packmol_unimplemented_keywords = Dict{String,String}(
)

# Global keywords (see `_set_input_keyword!` and `_read_packmol_input`)
const _global_keywords = (
    "filetype", "output", "tolerance", "radscale", "radscale_decay_loops", "writeout", "discale", "tolerance_precision",
    "constraint_precision", "maxit", "nloop", "movefrac", "nthreads", "periodic_boundary_style", "max_random_init",
    "adjust_constraints_on_init", "add_amber_ter", "amber_ter_preserve", "add_box_sides",
    "ignore_conect", "seed", "avoid_overlap", "screen_log",
    "check_initial_placement", "check", "continue_failed_initial_packing", "pbc", "unitcell",
    "voxel_place",
)

# Keywords of a `structure ... end structure` block with no use in Packmol.jl
# (accepted with a warning, as `packmol_legacy_keywords`).
const packmol_legacy_structure_keywords = Dict{String,String}(
    "nloop" => "nloop inside a structure block is a deprecated keyword and was ignored (use the global nloop).",
    "nloop0" => "nloop0 is a deprecated keyword and was ignored.",
    "restart_to" => "restart_to is a deprecated keyword and was ignored: to restart, give the output file (PDB, XYZ or mmCIF) to restart_from.",
)

# Input file line without its comment (anything after a `#`, as in Fortran
# Packmol) and surrounding whitespace.
_strip_comment(line::AbstractString) = strip(first(split(line, '#'; limit=2)))

"""
    read_packmol_input(input_file::String)
    read_packmol_input(io::IO)

Reads a Packmol input file and returns the `PackmolSystem` it describes, not packed
(pack it with `packmol(system)`). Paths in the input are relative to the directory of
`input_file`, or to the current directory if the input is read from an `IO` stream
(e.g. `stdin`).
"""
read_packmol_input(input::Union{String,IO}; D::Int=3, T::DataType=Float64) =
    read_packmol_input(input, Val(D), T)

function read_packmol_input(input_file::String, ::Val{D}, ::Type{T}) where {D,T}
    input_file = expanduser(input_file)
    return _read_packmol_input(open(_read_lines, input_file), input_file, dirname(input_file), Val(D), T)
end

read_packmol_input(io::IO, ::Val{D}, ::Type{T}) where {D,T} =
    _read_packmol_input(_read_lines(io), "", "", Val(D), T)

# `readlines` is not compatible with `juliac --trim`
function _read_lines(io::IO)
    lines = String[]
    while !eof(io)
        push!(lines, readline(io))
    end
    return lines
end

# With `screen_log no`, the warnings of the parser are silenced too, also those of the
# lines before the keyword
function _read_packmol_input(lines::Vector{String}, input_file::String, input_dir::String, ::Val{D}, ::Type{T}) where {D,T}
    screen_log = !any(_is_screen_log_no, lines)
    screen_log_before = _screen_log[]
    _screen_log[] = screen_log_before && screen_log
    try
        return _read_packmol_input_lines(lines, input_file, input_dir, Val(D), T)
    finally
        _screen_log[] = screen_log_before
    end
end
function _is_screen_log_no(line::AbstractString)
    words = split(_strip_comment(line))
    return length(words) == 2 && words[1] == "screen_log" && words[2] == "no"
end

function _read_packmol_input_lines(lines::Vector{String}, input_file::String, input_dir::String, ::Val{D}, ::Type{T}) where {D,T}
    sys = PackmolSystem{D,T}(; input_file, output_file="")
    structure_section = false
    for line in lines
        line = _strip_comment(line)
        isempty(line) && continue
        words = split(line)
        keyword, values = words[1], @view(words[2:end])
        # `structure`/`end structure` must be recognized regardless of
        # section state (to track it in the first place); everything
        # else inside a structure block is skipped here unconditionally
        # and left entirely to the second pass (`read_structure_data`
        # below) — otherwise a keyword name valid both globally and
        # per-structure (e.g. `restart_from`) would be wrongly applied
        # system-wide from inside a `structure ... end structure` block.
        if keyword == "structure"
            structure_section = true
            continue
        end
        if keyword == "end" && !isempty(values) && values[1] == "structure"
            structure_section = false
            continue
        end
        if structure_section
            continue
        end
        # Value-less flag keywords (e.g. `add_box_sides`)
        # appear on their own line with nothing after them; an empty string
        # is passed through for those.
        if _set_input_keyword!(sys, keyword, isempty(values) ? "" : first(values))
            continue
        elseif haskey(packmol_legacy_keywords, keyword)
            _warn(packmol_legacy_keywords[keyword])
            continue
        elseif haskey(packmol_unimplemented_keywords, keyword)
            _warn(packmol_unimplemented_keywords[keyword])
            continue
        end
        if keyword == "pbc"
            _set_pbc!(sys, values)
            continue
        end
        if keyword == "unitcell"
            _set_unitcell!(sys, values)
            continue
        end
        throw(ArgumentError("Keyword $keyword not recognized"))
    end
    if structure_section
        throw(ArgumentError("Structure block not closed with 'end structure'"))
    end
    # As the structure files, the other files named in the input file are
    # relative to its directory (an absolute path is kept as is by `joinpath`).
    isempty(sys.output_file) || (sys.output_file = joinpath(input_dir, sys.output_file))
    restart_from = sys.restart_from
    restart_from isa String && (sys.restart_from = joinpath(input_dir, restart_from))
    # Read structure data
    structure_input = IOBuffer()
    for line in lines
        line = _strip_comment(line)
        isempty(line) && continue
        words = split(line)
        if words[1] == "structure"
            structure_section = true
            structure_input = IOBuffer()
            print(structure_input, line, "\n")
            continue
        end
        if words[1] == "end" && length(words) > 1 && words[2] == "structure"
            print(structure_input, line, "\n")
            append!(sys.structure_types,
                _read_structure_types(structure_input, sys.tolerance, Val(D), T; path=input_dir, filetype=sys.filetype)
            )
            continue
        end
        if structure_section
            print(structure_input, line, "\n")
            continue
        end
    end
    unitcell, unitcell_center = sys.unitcell, sys.unitcell_center
    if !isnothing(unitcell) && !isnothing(unitcell_center)
        sys.wigner_seitz_cell = _is_wigner_seitz_cell(unitcell)
    end
    _check_fixed_inside_cell(sys)
    _check_restart_files(sys)
    sys.periodic_boundary_style != :triclinic && isnothing(sys.unitcell) && throw(ArgumentError("""\n
        periodic_boundary_style requires periodic boundary conditions (pbc or unitcell).

    """))
    #
    # Initialize atom data and molecule position arrays
    #
    sys.atoms, sys.molecule_positions, sys.nmols = _atoms_and_molecule_positions(sys.structure_types)
    return sys
end
