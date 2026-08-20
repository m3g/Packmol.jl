
@kwdef mutable struct PackmolSystem{D,T}
    filetype::String = "pdb"
    input_file::String
    output_file::String
    tolerance::T = 2.0
    # Loose-start radius scale factor for the packing objective (Fortran Packmol's
    # `discale`): atom radii are inflated by this factor early in packing, then
    # decayed toward 1.0 as improvement stalls (see packmol_main.jl).
    radscale::T = 1.2
    structure_types::Vector{StructureType{D,T}} = StructureType{D,T}[]
    tolerance_precision::T = 1e-2
    constraint_precision::T = 1e-2
    # Relative-improvement threshold below which the stall detector counts a
    # plateau iteration in either the minimum interatomic distance or the
    # maximum constraint violation during a packing chunk (see
    # packmol_main.jl).
    stall_tolerance::T = 1e-2 * tolerance_precision
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
    restart_from::Union{Nothing,String,Vector{<:Atom}} = nothing
    restart_to::Union{Nothing,String} = nothing
    # Unit cell for periodic boundary conditions (nothing = no PBC)
    unitcell::Union{Nothing, Matrix{T}} = nothing
    # Reference center for PBC wrapping (constraints evaluated relative to this point)
    unitcell_center::Union{Nothing, SVector{D,T}} = nothing
    # Internal data for the optimization
    nmols::Int = 0
    atoms::Vector{AtomData{T}} = AtomData{T}[]
    molecule_positions::Vector{MoleculePosition{D,T}} = MoleculePosition{D,T}[]
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
    atoms, molecule_positions, nmols = _atoms_and_molecule_positions(structure_types)
    return PackmolSystem{D,T}(;
        structure_types, atoms, molecule_positions, nmols,
        input_file, output_file=output, kargs...,
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
function _parse_value(T::DataType, keyword::String, input_value; _val_check=x -> nothing)
    input_value isa T && return input_value
    (T == String && input_value isa AbstractString) && return string(input_value)
    value = tryparse(T, input_value)
    if isnothing(value)
        throw(ArgumentError("""\n
            Could not parse value $input_value for keyword $keyword. Expected type: $T.

        """))
    end
    _val_check(value)
    return value
end

function _parse_options(T, name::String, val, options::Tuple)
    for pair in options
        if first(pair) == val
            return last(pair)
        end
    end
    throw(ArgumentError("""\n
        Value "$val" is not valid for option: $name.

    """))
end

#! format: off
packmol_input_keywords = Dict{String,Function}(
    "filetype"                 => (T, val) -> (:filetype, _parse_value(String, "filetype", val)),
    "output"                   => (T, val) -> (:output_file, _parse_value(String, "output", val)),
    "tolerance"                => (T, val) -> (:tolerance, _parse_value(T, "tolerance", val)),
    "radscale"                 => (T, val) -> (:radscale, _parse_value(T, "radscale", val)),
    "discale"                  => (T, val) -> begin
        @warn "discale is a legacy keyword name; use radscale instead." maxlog=1
        (:radscale, _parse_value(T, "discale", val))
    end,
    "tolerance_precision"      => (T, val) -> (:tolerance_precision, _parse_value(T, "tolerance_precision", val)),
    "constraint_precision"     => (T, val) -> (:constraint_precision, _parse_value(T, "constraint_precision", val)),
    "stall_tolerance"          => (T, val) -> (:stall_tolerance, _parse_value(T, "stall_tolerance", val)),
    "maxit"                    => (T, val) -> (:maxit, _parse_value(Int, "max_iter", val)),
    "max_random_init"          => (T, val) -> (:max_random_init, _parse_value(Int, "max_random_init", val)),
    "adjust_constraints_on_init" => (T, val) -> (:adjust_constraints_on_init, _parse_options(String, "adjust_constraints_on_init", val, ("yes" => true, "no" => false))),
    "add_amber_ter"            => (T, val) -> (:add_amber_ter, _parse_value(Bool, "add_amber_ter", val)),
    "amber_ter_preserve"       => (T, val) -> (:amber_ter_preserve, true),
    "add_box_sides"            => (T, val) -> (:add_box_sides, true),
    "connect"                  => (T, val) -> (:connect, _parse_options(String, "connect", val, ("yes" => true, "no" => false))),
    "randominitialpoint"       => (T, val) -> (:random_initial_point, true),
    "seed"                     => (T, val) -> (:seed, _parse_value(Int, "seed", val)),
    "avoid_overlap"            => (T, val) -> (:avoid_overlap, _parse_options(String, "avoid_overlap", val, ("yes" => true, "no" => false))),
    "writeout"                 => (T, val) -> (:writeout, _parse_value(Int, "writeout", val)),
    "writebad"                 => (T, val) -> (:writebad, true),
    "optimization_print_level" => (T, val) -> (:optimization_print_level, _parse_value(Int, "optim_print_level", val)),
    "chkgrad"                  => (T, val) -> (:chkgrad, true),
    "restart_from"             => (T, val) -> (:restart_from, _parse_value(String, "restart_from", val)),
    "restart_to"               => (T, val) -> (:restart_to, _parse_value(String, "restart_to", val)),
)
#! format: on

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

packmol_legacy_keywords = Dict{String,String}(
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
# here, rather than in packmol_input_keywords, so that using them produces a
# clear warning instead of a crash while they remain unimplemented.
packmol_unimplemented_keywords = Dict{String,String}(
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
function read_packmol_input(input_file::String; D::Int=3, T::DataType=Float64)
    input_data = Dict{Symbol,Any}(
        :input_file => input_file,
        :structure_types => StructureType{D,T}[],
    )
    structure_section = false
    open(input_file) do io
        for line in eachline(io)
            line = strip(line)
            (startswith(line, "#") || isempty(line)) && continue
            keyword, values... = split(line)
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
            # Handle no-value keywords first
            if keyword == "check"
                input_data[:check] = true
                continue
            end
            if haskey(packmol_input_keywords, keyword)
                # Value-less flag keywords (e.g. `randominitialpoint`, `add_box_sides`)
                # appear on their own line with nothing after them; their closures
                # ignore `val`, so an empty string is passed through for those.
                val = isempty(values) ? "" : first(values)
                field_name, field_value = packmol_input_keywords[keyword](T, val)
                input_data[field_name] = field_value
                continue
            elseif haskey(packmol_legacy_keywords, keyword)
                @warn begin
                    packmol_legacy_keywords[keyword]
                end _line = line _file = input_file
                continue
            elseif haskey(packmol_unimplemented_keywords, keyword)
                @warn begin
                    packmol_unimplemented_keywords[keyword]
                end _line = line _file = input_file
                continue
            end
            # pbc: orthorhombic PBC keyword
            #   3 values (side lengths): center at origin
            #   6 values (xmin ymin zmin xmax ymax zmax): center at midpoint
            if keyword == "pbc"
                vals = [_parse_value(T, "pbc", v) for v in values]
                if length(vals) == D
                    # 3 values: side lengths, center at origin
                    input_data[:unitcell] = Matrix{T}(Diagonal(vals))
                    input_data[:unitcell_center] = zero(SVector{D,T})
                elseif length(vals) == 2 * D
                    # 6 values: xmin ymin zmin xmax ymax zmax
                    lo = SVector{D,T}(vals[1:D]...)
                    hi = SVector{D,T}(vals[D+1:2*D]...)
                    sides = hi - lo
                    input_data[:unitcell] = Matrix{T}(Diagonal(Vector(sides)))
                    input_data[:unitcell_center] = (lo + hi) / 2
                else
                    throw(ArgumentError("pbc keyword requires $D (side lengths) or $(2*D) (min/max coords) values, got $(length(vals))"))
                end
                continue
            end
            # unitcell: general PBC keyword (6 values: a b c α β γ in CRYST1 convention)
            # Center at origin.
            if keyword == "unitcell"
                if length(values) != 2 * D
                    throw(ArgumentError("unitcell keyword requires $(2*D) values (a b [c] α β [γ]), got $(length(values))"))
                end
                params = [_parse_value(T, "unitcell", v) for v in values]
                input_data[:unitcell] = Matrix{T}(unitcell_matrix(T, params...))
                input_data[:unitcell_center] = zero(SVector{D,T})
                continue
            end
            throw(ArgumentError("Keyword $keyword not recognized"))
        end
    end
    if structure_section
        throw(ArgumentError("Structure block not closed with 'end structure'"))
    end
    # Read structure data
    open(input_file) do io
        local structure_input
        for line in eachline(io)
            line = strip(line)
            (startswith(line, "#") || isempty(line)) && continue
            keyword, values... = split(line)
            if keyword == "structure"
                structure_section = true
                structure_input = IOBuffer()
                print(structure_input, line, "\n")
                continue
            end
            if keyword == "end" && values[1] == "structure"
                print(structure_input, line, "\n")
                push!(input_data[:structure_types], 
                    read_structure_data(structure_input, input_data[:tolerance]; D, T, path=dirname(input_file))
                )
                continue
            end
            if structure_section
                print(structure_input, line, "\n")
                continue
            end
        end
    end
    #
    # Initialize atom data and molecule position arrays
    #
    atoms, molecule_positions, nmols = _atoms_and_molecule_positions(input_data[:structure_types])
    input_data[:atoms] = atoms
    input_data[:molecule_positions] = molecule_positions
    input_data[:nmols] = nmols

    return PackmolSystem{D,T}(; input_data...)
end
