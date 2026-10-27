
@kwdef mutable struct PackmolSystem{D,T}
    filetype::String = "pdb"
    input_file::String
    output_file::String
    tolerance::T = 2.0
    # Loose-start radius scale factor for the packing objective (Fortran Packmol's
    # `discale`): atom radii are inflated by this factor at the start of packing
    # (loop 0), then decayed exponentially toward 1.0 by loop `nloop`, the
    # same for every atom regardless of movebad! relocations (see
    # `_radscale_at` and the outer packing loop in packmol_main.jl).
    radscale::T = 1.2
    structure_types::Vector{StructureType{D,T}} = StructureType{D,T}[]
    tolerance_precision::T = 1e-1
    constraint_precision::T = 1e-2
    # Relative-improvement threshold below which the stall detector counts a
    # plateau iteration in either the minimum interatomic distance or the
    # maximum constraint violation during a packing chunk (see
    # packmol_main.jl).
    stall_tolerance::T = 1e-3 * tolerance_precision
    max_iter::Int = 1000
    max_random_init::Int=20
    adjust_constraints_on_init::Bool=true
    add_amber_ter::Bool = false
    amber_ter_preserve::Bool = false
    add_box_sides::Bool = false
    connect::Bool = false
    random_initial_point::Bool = false
    seed::Int = 1234567
    avoid_overlap::Bool = true
    writeout::Int = 20
    writebad::Bool = false
    optim_print_level::Int = 0
    chkgrad::Bool = false
    check::Bool = false
    # Whole-system restart (Fortran's `restart_from`/`restart_to` with no
    # structure type, i.e. index 0): `nothing` means unset. When
    # `restart_from` is set, packmol() skips the whole initial-approximation
    # pipeline entirely and reads every molecule's position from this
    # source instead — the point of a whole-system restart is to skip that
    # (potentially expensive, for a large system) phase, not just to get
    # the same final starting point some other way.
    # `restart_from` accepts either a file path (a raw Packmol-format
    # restart file, or a `.pdb` — recognized by extension — whose atoms are
    # rigid-body-aligned against each structure type's own template to
    # recover (cm, angles); see `_align_molecule` in write_output.jl) or,
    # Julia-API only, a `Vector{<:Atom}` already in memory (skipping the
    # file entirely). `restart_to` only ever writes the raw format.
    # See also `StructureType.restart_from`/`restart_to` for the
    # per-structure-type equivalents.
    restart_from::Union{Nothing,String,Vector{Atom{Nothing}}} = nothing
    restart_to::Union{Nothing,String} = nothing
    # Unit cell for periodic boundary conditions (nothing = no PBC)
    unitcell::Union{Nothing, Matrix{T}} = nothing
    # Reference center for PBC wrapping (constraints evaluated relative to this point)
    unitcell_center::Union{Nothing, SVector{D,T}} = nothing
    # Which fundamental-domain shape `get_atoms`/`write_output` wrap molecule
    # positions into, when `unitcell` is set: `:triclinic` (the ordinary
    # skewed-parallelepiped shape, e.g. what CRYST1/PDB output implies),
    # `:dodecahedral` (the rhombic-dodecahedron/Wigner-Seitz "compact"
    # shape), or `:octahedral` (the truncated-octahedron/Wigner-Seitz
    # "compact" shape). Set via `triclinic_to_dodecahedral(::PackmolSystem)`/
    # `dodecahedral_to_triclinic(::PackmolSystem)` (or their `octahedral`
    # counterparts) rather than directly.
    periodic_boundary_style::Symbol = :triclinic
    # Internal data for the optimization
    nmols::Int = 0
    atoms::Vector{AtomData{T}} = AtomData{T}[]
    molecule_positions::Vector{MoleculePosition{D,T}} = MoleculePosition{D,T}[]
    # Outcome of the last `packmol(::PackmolSystem)` call on this system:
    # `:not_packed` (never packed, e.g. freshly built, or only `check`-mode
    # placement), `:packing_ready` (converged within `tolerance_precision`/
    # `constraint_precision`), or `:failed_packing` (ran out of `nloop`
    # loops short of that). `packmol` returns this same `PackmolSystem`
    # (mutated in place) rather than a bare `Bool`, precisely so this field
    # is how callers check the outcome — see `packmol(::PackmolSystem)`.
    status::Symbol = :not_packed
end

#
# Build the per-atom (`AtomData`) and per-molecule (`MoleculePosition`) arrays
# from a vector of `StructureType`s: one is a flat list of every atom across
# every molecule of every structure type, tagging it with its molecule index,
# structure type index, radius, and constraint indices; the other is a
# placeholder (zeroed) position per molecule, later filled in by
# `initialize_molecules!`. Shared by `read_packmol_input` and the
# `PackmolSystem` outer constructor below.
#
function _atoms_and_molecule_positions(structure_types::Vector{StructureType{D,T}}) where {D,T}
    mol_index = 0
    atoms = AtomData{T}[]
    molecule_positions = MoleculePosition{D,T}[]
    for (itype, st) in enumerate(structure_types)
        for _ in 1:st.number_of_molecules
            mol_index += 1
            push!(molecule_positions, MoleculePosition(zeros(SVector{D,T}), zeros(SVector{D,T})))
            for iatom in 1:st.natoms
                push!(atoms, AtomData(mol_index, itype, st.radii[iatom], st.atom_constraints[iatom]))
            end
        end
    end
    return atoms, molecule_positions, mol_index
end

#
# Append one constraint to a structure type's own constraint list, wiring it
# into every atom's `atom_constraints` (not just the atoms present when the
# structure type was originally built) — including atoms whose own
# `atoms ... end atoms` block already narrowed them to a custom subset of
# constraints, since this one is a whole-system requirement, not a per-atom
# feature toggle.
#
function _push_constraint!(st::StructureType{D,T}, c::AnyConstraint{T}) where {D,T}
    push!(st.constraints, c)
    idx = length(st.constraints)
    for ac in st.atom_constraints
        push!(ac, idx)
    end
    return st
end

#
# Implicitly confine every free (non-fixed) structure type to (a slightly
# inflated version of) the PBC unit cell, whenever PBC is active.
#
# Without this, PBC only affects interatomic distances (via CellListMap) and
# how *explicit* constraints get evaluated (wrapped into one cell image) —
# nothing stops a molecule's own (unbounded) center-of-mass optimization
# variable from drifting freely along any direction no explicit constraint
# happens to bound. That's mostly harmless for a structure type with its own
# `inside box`/`inside cube` already matching the cell, but a constraint like
# `above plane`/`below plane` is a half-space with no periodicity of its
# own — nothing else ties it to one well-defined periodic image, and without
# some other bound the molecule is free to wander along the plane's
# unconstrained directions (or its own unconstrained side).
#
# For an orthorhombic cell (diagonal `unitcell`) this reduces to a single
# `InsideBox` spanning the cell; the general (triclinic) case needs one
# inward-facing `Plane` pair per pair of opposing cell faces instead, since a
# box's axis-aligned edges can't represent a sheared cell.
#
# Each face is pushed outward by `tolerance / 2` (so the box/planes end up
# `tolerance` larger than the literal cell, in total, along each axis): an
# atom sitting exactly on the true cell boundary needs its own radius'
# worth of room (the default atom radius is `tolerance / 2`) to interact
# naturally with its periodic image across that boundary, the same way it
# would with a same-radius neighbor anywhere in the interior. Without this
# slack, the implicit confinement would push back *before* that natural
# interaction has a chance to happen, visibly flattening the packing right
# at the cell boundary instead of letting it look like an unbroken,
# continuous interior.
function _add_pbc_confinement_constraints!(
    structure_types::Vector{StructureType{D,T}},
    unitcell::Matrix{T},
    unitcell_center::SVector{D,T},
    tolerance::T,
) where {D,T}
    is_orthorhombic = D == 3 && all(i == j || unitcell[i, j] == zero(T) for i in 1:D, j in 1:D)
    if is_orthorhombic
        sides = SVector{D,T}(ntuple(i -> unitcell[i, i] + tolerance, D))
        confinement = (InsideBox(unitcell_center, sides),)
    else
        a = SVector{3,T}(unitcell[1, 1], unitcell[2, 1], unitcell[3, 1])
        b = SVector{3,T}(unitcell[1, 2], unitcell[2, 2], unitcell[3, 2])
        c = SVector{3,T}(unitcell[1, 3], unitcell[2, 3], unitcell[3, 3])
        confinement = ntuple(3) do k
            u, v, w = (a, b, c)[k], (b, c, a)[k], (c, a, b)[k]
            n = cross(v, w)
            n = n / norm(n)
            dot(u, n) < zero(T) && (n = -n)
            d0 = dot(n, unitcell_center)
            half = dot(n, u) / 2 + tolerance / 2
            (BelowPlane(n, d0 + half), AbovePlane(n, d0 - half))
        end
        confinement = (confinement[1]..., confinement[2]..., confinement[3]...)
    end
    for st in structure_types
        st.fixed.fixed && continue
        for c in confinement
            _push_constraint!(st, c)
        end
    end
    return structure_types
end

"""
    PackmolSystem(structure_types::Vector{<:StructureType}; output::String, kargs...)

Build a `PackmolSystem` directly from a vector of `StructureType`s (see
[`structure_type`](@ref)), without going through an input file. `output` sets
the output file path; any other `PackmolSystem` field (`tolerance`, `seed`,
`unitcell`, ...) can be given as a keyword argument.
"""
function PackmolSystem(
    structure_types::Vector{StructureType{D,T}};
    output::String,
    input_file::String = "",
    kargs...,
) where {D,T}
    # Implicitly confine every free structure type to (a slightly inflated
    # version of) the PBC cell, when PBC is active — see
    # `_add_pbc_confinement_constraints!`'s own docstring for why.
    unitcell_raw = get(kargs, :unitcell, nothing)
    if !isnothing(unitcell_raw)
        unitcell = Matrix{T}(unitcell_raw)
        unitcell_center = SVector{D,T}(get(kargs, :unitcell_center, zero(SVector{D,T})))
        tolerance = T(get(kargs, :tolerance, T(2.0)))
        _add_pbc_confinement_constraints!(structure_types, unitcell, unitcell_center, tolerance)
    end
    atoms, molecule_positions, nmols = _atoms_and_molecule_positions(structure_types)
    kw = values(kargs)
    if haskey(kw, :restart_from)
        kw = merge(kw, (; restart_from=_restart_source(kw.restart_from)))
    end
    return PackmolSystem{D,T}(;
        structure_types, atoms, molecule_positions, nmols,
        input_file, output_file=output, kw...,
    )
end

function _indent(s::AbstractString; n=4)
    indented_str = IOBuffer()
    for line in eachline(IOBuffer(s))
        println(indented_str, repeat(" ", n) * line)
    end
    return String(take!(indented_str))
end

function Base.show(io::IO, ::MIME"text/plain", sys::PackmolSystem{D,T}) where {D,T}
    printstyled(io, "PackmolSystem{$D,$T}"; bold=true, color=:blue)
    if length(sys.input_file) > 0 
        printstyled(" - read from: $(basename(sys.input_file))"; color=:blue) 
    end
    println(io)
    printstyled(io, _indent("Number of types of structures: $(length(sys.structure_types))\n"); bold=true)
    printstyled(io, _indent("Total number of molecules: $(sys.nmols)\n"); bold=true)
    printstyled(io, _indent("Structure types:\n"); bold=true)
    for (i,st) in enumerate(sys.structure_types)
        print(io, _indent("$i."*_show(st)))
    end
    printstyled(io, _indent("Options:\n"); bold=true)
    for field in fieldnames(PackmolSystem)
        if !(field in (:structure_types, :input_file, :atoms, :molecule_positions, :nmols))
            print(io, _indent("$field: $(getfield(sys, field))"; n=8))
        end
    end
end

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
_warn(msg::AbstractString) = _with_stderr(io -> println(io, "Warning: ", msg))

#
# Set the field of `sys` corresponding to a global keyword of the input file,
# with value `val` (an empty string for the value-less flag keywords, e.g.
# `randominitialpoint`). Returns `false` if `keyword` is not one of them. An
# explicit chain of comparisons, rather than a dictionary of closures, so
# that each field is set with a value of its own (inferred) type.
#
function _set_input_keyword!(sys::PackmolSystem{D,T}, keyword::AbstractString, val::AbstractString) where {D,T}
    if keyword == "filetype"
        sys.filetype = _parse_value(String, "filetype", val)
    elseif keyword == "output"
        sys.output_file = expanduser(_parse_value(String, "output", val))
    elseif keyword == "tolerance"
        sys.tolerance = _parse_value(T, "tolerance", val)
    elseif keyword == "radscale"
        sys.radscale = _parse_value(T, "radscale", val)
    elseif keyword == "discale"
        _warn("discale is a legacy keyword name; use radscale instead.")
        sys.radscale = _parse_value(T, "discale", val)
    elseif keyword == "tolerance_precision"
        sys.tolerance_precision = _parse_value(T, "tolerance_precision", val)
    elseif keyword == "constraint_precision"
        sys.constraint_precision = _parse_value(T, "constraint_precision", val)
    elseif keyword == "stall_tolerance"
        sys.stall_tolerance = _parse_value(T, "stall_tolerance", val)
    elseif keyword == "maxit"
        sys.max_iter = _parse_value(Int, "maxit", val)
    elseif keyword == "max_random_init"
        sys.max_random_init = _parse_value(Int, "max_random_init", val)
    elseif keyword == "adjust_constraints_on_init"
        sys.adjust_constraints_on_init = _parse_options("adjust_constraints_on_init", val, _yes_no)
    elseif keyword == "add_amber_ter"
        sys.add_amber_ter = true
    elseif keyword == "amber_ter_preserve"
        sys.amber_ter_preserve = true
    elseif keyword == "add_box_sides"
        sys.add_box_sides = true
    elseif keyword == "connect"
        sys.connect = _parse_options("connect", val, _yes_no)
    elseif keyword == "randominitialpoint"
        sys.random_initial_point = true
    elseif keyword == "seed"
        sys.seed = _parse_value(Int, "seed", val)
    elseif keyword == "avoid_overlap"
        sys.avoid_overlap = _parse_options("avoid_overlap", val, _yes_no)
    elseif keyword == "writeout"
        sys.writeout = _parse_value(Int, "writeout", val)
    elseif keyword == "writebad"
        sys.writebad = true
    elseif keyword == "optimization_print_level"
        sys.optim_print_level = _parse_value(Int, "optimization_print_level", val)
    elseif keyword == "chkgrad"
        sys.chkgrad = true
    elseif keyword == "check"
        sys.check = true
    elseif keyword == "restart_from"
        sys.restart_from = expanduser(_parse_value(String, "restart_from", val))
    elseif keyword == "restart_to"
        sys.restart_to = expanduser(_parse_value(String, "restart_to", val))
    else
        return false
    end
    return true
end

#=
    unitcell_matrix(a, b, c, α, β, γ)

Convert CRYST1-style unit cell parameters (sides a, b, c and angles α, β, γ in degrees)
to a 3×3 unit cell matrix (columns are the cell vectors). Uses the PDB convention:
  a along x, b in the xy-plane.
=#
function unitcell_matrix(::Type{T}, a, b, c, α_deg, β_deg, γ_deg) where {T}
    α = T(α_deg * π / 180)
    β = T(β_deg * π / 180)
    γ = T(γ_deg * π / 180)
    ax = T(a)
    bx = T(b) * cos(γ)
    by = T(b) * sin(γ)
    cx = T(c) * cos(β)
    cy = T(c) * (cos(α) - cos(β) * cos(γ)) / sin(γ)
    cz = sqrt(T(c)^2 - cx^2 - cy^2)
    @SMatrix [
        ax  bx  cx
        zero(T)  by  cy
        zero(T) zero(T) cz
    ]
end

const packmol_legacy_keywords = Dict{String,String}(
    "fscale" => "fscale legacy keyword was ignored.",
    "fbins" => "fbins legacy keyword was ignored.",
    "iprint1" => "iprint1 legacy keyword was ignored, instead use: optim_print_level",
    "iprint2" => "iprint1 legacy keyword was ignored, instead use: optim_print_level",
    "iprint3" => "iprint1 legacy keyword was ignored, instead use: optim_print_level",
    "precision" => "precision legacy keyword was ignored, instead use: tolerance_precision and/or constraint_precision",
    "packall" => "packall is ignored and the only option",
)

# Keywords recognized by Fortran Packmol but not yet wired into Packmol.jl
# (see the Phase 1 "Missing keyword parsing" checklist in CLAUDE.md). Listed
# here, rather than in `_set_input_keyword!`, so that using them produces a
# clear warning instead of a crash while they remain unimplemented.
const packmol_unimplemented_keywords = Dict{String,String}(
    "movefrac" => "movefrac keyword is not yet implemented in Packmol.jl and was ignored.",
    "movebadrandom" => "movebadrandom keyword is not yet implemented in Packmol.jl and was ignored.",
    "use_short_tol" => "use_short_tol keyword is not yet implemented in Packmol.jl and was ignored.",
    "short_tol_dist" => "short_tol_dist keyword is not yet implemented in Packmol.jl and was ignored.",
    "short_tol_scale" => "short_tol_scale keyword is not yet implemented in Packmol.jl and was ignored.",
    "short_radius" => "short_radius keyword is not yet implemented in Packmol.jl and was ignored.",
    "short_radius_scale" => "short_radius_scale keyword is not yet implemented in Packmol.jl and was ignored.",
    "nloop" => "nloop keyword is not yet implemented in Packmol.jl and was ignored.",
)

#=

    read_packmol_input

Reads the packmol input file to create a `PackmolSystem` object.

=#
read_packmol_input(input_file::String; D::Int=3, T::DataType=Float64) =
    read_packmol_input(input_file, Val(D), T)

function read_packmol_input(input_file::String, ::Val{D}, ::Type{T}) where {D,T}
    input_file = expanduser(input_file)
    sys = PackmolSystem{D,T}(; input_file, output_file="")
    stall_tolerance_set = false
    structure_section = false
    for line in eachline(input_file)
        line = strip(line)
        (startswith(line, "#") || isempty(line)) && continue
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
        # Value-less flag keywords (e.g. `randominitialpoint`, `add_box_sides`)
        # appear on their own line with nothing after them; an empty string
        # is passed through for those.
        if _set_input_keyword!(sys, keyword, isempty(values) ? "" : first(values))
            stall_tolerance_set |= keyword == "stall_tolerance"
            continue
        elseif haskey(packmol_legacy_keywords, keyword)
            _warn(packmol_legacy_keywords[keyword])
            continue
        elseif haskey(packmol_unimplemented_keywords, keyword)
            _warn(packmol_unimplemented_keywords[keyword])
            continue
        end
        # pbc: orthorhombic PBC keyword
        #   3 values (side lengths): center at origin
        #   6 values (xmin ymin zmin xmax ymax zmax): center at midpoint
        if keyword == "pbc"
            if !isempty(values) && values[1] == "dodecahedral"
                unitcell, unitcell_center = parse_pbc_dodecahedral(T, values[2:end], Val(D))
                sys.unitcell, sys.unitcell_center = unitcell, unitcell_center
                continue
            end
            if !isempty(values) && values[1] == "octahedral"
                unitcell, unitcell_center = parse_pbc_octahedral(T, values[2:end], Val(D))
                sys.unitcell, sys.unitcell_center = unitcell, unitcell_center
                continue
            end
            vals = [_parse_value(T, "pbc", v) for v in values]
            if length(vals) == D
                # 3 values: side lengths, center at origin
                sys.unitcell = Matrix{T}(Diagonal(vals))
                sys.unitcell_center = zero(SVector{D,T})
            elseif length(vals) == 2 * D
                # 6 values: xmin ymin zmin xmax ymax zmax
                lo = SVector{D,T}(@view(vals[1:D]))
                hi = SVector{D,T}(@view(vals[D+1:2*D]))
                sides = hi - lo
                sys.unitcell = Matrix{T}(Diagonal(Vector(sides)))
                sys.unitcell_center = (lo + hi) / 2
            else
                throw(ArgumentError("pbc keyword requires $D (side lengths) or $(2*D) (min/max coords) values, got $(length(vals))"))
            end
            continue
        end
        # unitcell: general PBC keyword (6 values: a b c α β γ in CRYST1 convention)
        # Center at origin.
        if keyword == "unitcell"
            if D != 3 || length(values) != 6
                throw(ArgumentError("unitcell keyword requires 6 values (a b c α β γ), got $(length(values))"))
            end
            p = [_parse_value(T, "unitcell", v) for v in values]
            sys.unitcell = Matrix{T}(unitcell_matrix(T, p[1], p[2], p[3], p[4], p[5], p[6]))
            sys.unitcell_center = zero(SVector{D,T})
            continue
        end
        throw(ArgumentError("Keyword $keyword not recognized"))
    end
    if structure_section
        throw(ArgumentError("Structure block not closed with 'end structure'"))
    end
    # The default stall tolerance is relative to the tolerance precision, which
    # might have been set in the input file.
    if !stall_tolerance_set
        sys.stall_tolerance = T(1e-3) * sys.tolerance_precision
    end
    # As the structure files, the other files named in the input file are
    # relative to its directory (an absolute path is kept as is by `joinpath`).
    input_dir = dirname(input_file)
    isempty(sys.output_file) || (sys.output_file = joinpath(input_dir, sys.output_file))
    restart_from = sys.restart_from
    restart_from isa String && (sys.restart_from = joinpath(input_dir, restart_from))
    restart_to = sys.restart_to
    restart_to isa String && (sys.restart_to = joinpath(input_dir, restart_to))
    # Read structure data
    structure_input = IOBuffer()
    for line in eachline(input_file)
        line = strip(line)
        (startswith(line, "#") || isempty(line)) && continue
        words = split(line)
        if words[1] == "structure"
            structure_section = true
            structure_input = IOBuffer()
            print(structure_input, line, "\n")
            continue
        end
        if words[1] == "end" && length(words) > 1 && words[2] == "structure"
            print(structure_input, line, "\n")
            push!(sys.structure_types, 
                read_structure_data(structure_input, sys.tolerance, Val(D), T; path=input_dir)
            )
            continue
        end
        if structure_section
            print(structure_input, line, "\n")
            continue
        end
    end
    # Implicitly confine every free structure type to (a slightly inflated
    # version of) the PBC cell, when PBC is active — see
    # `_add_pbc_confinement_constraints!`'s own docstring for why.
    unitcell, unitcell_center = sys.unitcell, sys.unitcell_center
    if !isnothing(unitcell) && !isnothing(unitcell_center)
        _add_pbc_confinement_constraints!(sys.structure_types, unitcell, unitcell_center, sys.tolerance)
    end
    #
    # Initialize atom data and molecule position arrays
    #
    sys.atoms, sys.molecule_positions, sys.nmols = _atoms_and_molecule_positions(sys.structure_types)
    return sys
end
