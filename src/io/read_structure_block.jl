#
# The contents of a structure block of an input file (`structure ... end structure`): the
# structure file, the number of molecules, the constraints, and the other keywords.
#
read_structure_data(input_file_block::AbstractString, args...; kwargs...) = 
    read_structure_data(IOBuffer(input_file_block), args...; kwargs...)

read_structure_data(input_file_block::IOBuffer, tolerance; D::Int=3, T::DataType=Float64, path::String="", filetype::String="pdb") =
    read_structure_data(input_file_block, tolerance, Val(D), T; path, filetype)

#
# Keywords of a structure block that are not used here: deprecated keywords and global
# keywords (which Fortran Packmol also accepts inside structure blocks) are ignored with
# a warning, and anything else is an error (e.g. a misspelled keyword).
function _other_structure_keyword(keyword::AbstractString, filename::AbstractString)
    if haskey(packmol_legacy_keywords, keyword)
        _warn(packmol_legacy_keywords[keyword])
    elseif keyword in _global_keywords
        _warn("$keyword is a global keyword, and was ignored inside the structure block of $(basename(filename)).")
    else
        throw(ArgumentError("""\n
            Keyword $keyword not recognized, in the structure block of $(basename(filename)).

        """))
    end
    return nothing
end

# The segment identifier of the PDB format has at most 4 characters
_check_segid(::Nothing) = nothing
function _check_segid(segid::String)
    if length(segid) > 4
        throw(ArgumentError("""\n
            segid must have at most 4 characters, got: $segid

        """))
    end
    return segid
end

#
# `constrain_rotation <axis> <center> <halfwidth>` lines of a structure block, optionally
# followed by the two sets of atoms that define the vector of the molecule:
#
#     constrain_rotation z 0. 15.
#       atoms 1 2 3
#       to 25 26
#     end constrain_rotation
#
# A `constrain_rotation` line starts such a block when it is followed by an `atoms` line
# and a `to` line (an `atoms` block of per-atom constraints has no `to`). Returns the
# block without these lines (read by `read_structure_data`), and the constraints.
#
function _read_rotation_constraints(input_file_block::IOBuffer, ::Val{D}, ::Type{T}) where {D,T}
    seekstart(input_file_block)
    lines = _read_lines(input_file_block)
    words = [split(_strip_comment(line)) for line in lines]
    significant = [i for i in eachindex(words) if !isempty(words[i])]
    rotation_constraints = RotationConstraint{D,T}[]
    skip = falses(length(lines))
    k = 1
    while k <= length(significant)
        i = significant[k]
        w = words[i]
        k += 1
        w[1] == "constrain_rotation" || continue
        length(w) >= 4 || throw(ArgumentError("""\n
            constrain_rotation requires an axis, a center and a half-width (degrees): $(strip(lines[i]))

        """))
        skip[i] = true
        from, to = Int[], Int[]
        if k + 1 <= length(significant) && words[significant[k]][1] == "atoms" && words[significant[k+1]][1] == "to"
            from = [parse(Int, v) for v in @view(words[significant[k]][2:end])]
            to = [parse(Int, v) for v in @view(words[significant[k+1]][2:end])]
            iend = k + 2 <= length(significant) ? significant[k+2] : 0
            (iend > 0 && words[iend][1] == "end" && length(words[iend]) >= 2 && words[iend][2] == "constrain_rotation") ||
                throw(ArgumentError("""\n
                    constrain_rotation: the `atoms` and `to` lines must be followed by `end constrain_rotation`.

                """))
            skip[significant[k]] = skip[significant[k+1]] = skip[iend] = true
            k += 3
        end
        axis = w[2] == "x" ? :x : w[2] == "y" ? :y : w[2] == "z" ? :z : :invalid
        push!(rotation_constraints, RotationConstraint{D,T}(axis, parse(T, w[3]), parse(T, w[4]); from, to))
    end
    filtered = IOBuffer()
    for i in eachindex(lines)
        skip[i] || println(filtered, lines[i])
    end
    return seekstart(filtered), rotation_constraints
end

function _check_rotation_constraints(rotation_constraints, fixed, natoms::Int)
    isempty(rotation_constraints) && return nothing
    fixed.fixed && throw(ArgumentError("""\n
        constrain_rotation cannot be used with a fixed molecule.

    """))
    for rc in rotation_constraints, i in Iterators.flatten((rc.from, rc.to))
        (1 <= i <= natoms) || throw(ArgumentError("""\n
            constrain_rotation: atom index $i out of range (the molecule has $natoms atoms).

        """))
    end
    return nothing
end

function read_structure_data(input_file_block::IOBuffer, tolerance, ::Val{D}, ::Type{T}; path::String="", filetype::String="pdb") where {D,T}
    types = _read_structure_types(input_file_block, tolerance, Val(D), T; path, filetype)
    return length(types) == 1 ? only(types) : types
end

#
# The structure types of a structure block of an input file: one, or one per
# conformation if the structure file has several models (see `_conformation_types`).
#
function _read_structure_types(input_file_block::IOBuffer, tolerance, ::Val{D}, ::Type{T}; path::String="", filetype::String="pdb") where {D,T}
    filename = ""
    conformations = Vector{Atom{Nothing}}[]
    conformation_files = String[]
    atoms = Atom{Nothing}[]
    connectivity = nothing
    connect = true
    number_of_molecules = 0
    reference_coordinates = SVector{D,T}[]
    radii = T[]
    atom_constraints = Vector{Int}[]
    constraints = AnyConstraint{T}[]
    fixed = zero(FixedMoleculeData{D,T})
    center = :none # :none, :geometric, or :mass
    residue_numbering = -1
    chain = nothing
    changechains = false
    segid = nothing
    radscale = nothing
    input_file_block, rotation_constraints = _read_rotation_constraints(input_file_block, Val(D), T)
    restart_from = nothing
    # Read basic structure data first
    atoms_block = false
    iconstraint = 0
    seekstart(input_file_block)
    for line in eachline(input_file_block)
        words = split(_strip_comment(line))
        isempty(words) && continue
        keyword, values = words[1], @view(words[2:end])
        if haskey(packmol_legacy_structure_keywords, keyword)
            _warn(packmol_legacy_structure_keywords[keyword])
            continue
        end
        if keyword == "atoms"
            atoms_block = true
        end
        if atoms_block && keyword == "end"
            atoms_block = false
        end
        if keyword == "structure"
            filename = joinpath(path, expanduser(String(values[1])))
            conformations, conformation_files, connectivity = _read_conformations(filename, filetype)
            atoms = conformations[1]
            reference_coordinates = [_atom_position(at, Val(D), T) for at in atoms]
            radii = fill(T(tolerance / 2), length(atoms))
            atom_constraints = [Int[] for _ in eachindex(atoms)]
        elseif keyword == "number"
            number_of_molecules = parse(Int, values[1]) 
        elseif keyword == "fixed"
            vals = [parse(T, v) for v in values]
            cm = SVector{D,T}(ntuple(i -> vals[i], D))
            angles = SVector{D,T}(ntuple(i -> vals[D+i], D))
            fixed = FixedMoleculeData(true, MoleculePosition(cm, angles))
        elseif keyword == "center"
            center = :geometric
        elseif keyword == "centerofmass"
            center = :mass
        elseif keyword == "resnumbers"
            residue_numbering = parse(Int, values[1])
        elseif keyword == "chain"
            chain = only(values[1])
        elseif keyword == "changechains"
            changechains = true
        elseif keyword == "connect"
            connect = _parse_options("connect", isempty(values) ? "" : values[1], _yes_no)
        elseif keyword == "segid"
            segid = _check_segid(String(values[1]))
        elseif keyword == "radscale"
            radscale = _parse_value(T, "radscale", values[1])
        elseif keyword == "restart_from"
            restart_from = joinpath(path, expanduser(String(values[1])))
        elseif keyword in constraint_placements
            iconstraint += 1
            push_parsed_constraint!(constraints, "$keyword $(values[1])", @view(values[2:end]), filename)
            if !atoms_block
                for iatom in eachindex(atoms)
                    push!(atom_constraints[iatom], iconstraint)
                end
            end
        elseif !atoms_block && keyword == "radius"
            radii .= parse(T, values[1])
        elseif keyword in ("atoms", "end", "radius")
            # `atoms` blocks are read below
        else
            _other_structure_keyword(keyword, filename)
        end
    end
    # If molecule is fixed, apply transformation to obtain the reference coordinates
    _apply_fixed_center!(reference_coordinates, atoms, fixed, center)
    # A fixed structure is a single molecule (`number` is optional, as in Fortran Packmol)
    fixed.fixed && number_of_molecules == 0 && (number_of_molecules = 1)
    #
    # Read custom atom radii and set constraints, according to atom blocks
    #
    atoms_block = false
    iconstraint = 0
    atoms_in_block = Int[]
    seekstart(input_file_block)
    for line in eachline(input_file_block)
        words = split(_strip_comment(line))
        isempty(words) && continue
        keyword, values = words[1], @view(words[2:end])
        if keyword == "atoms" 
            atoms_block = true
            atoms_in_block = [parse(Int, v) for v in values]
        end
        if keyword == "end" && atoms_block
            atoms_block = false
        end
        if atoms_block && keyword == "radius"
            radii[atoms_in_block] .= parse(T, values[1])
        end
        if keyword in constraint_placements
            iconstraint += 1
            if atoms_block
                for iatom in atoms_in_block
                    push!(atom_constraints[iatom], iconstraint)
                end
            end
        end
    end
    if changechains && !isnothing(chain)
        error("'changechains' and 'chain' input parameters are not compatible for a single structure.")
    end
    _check_rotation_constraints(rotation_constraints, fixed, length(atoms))
    st = StructureType{D,T}(;
        filename, natoms=length(atoms), atoms, number_of_molecules, fixed, reference_coordinates, radii,
        residue_numbering, chain, changechains, segid, radscale, rotation_constraints, restart_from,
        connectivity=connect ? connectivity : nothing,
        constraints, atom_constraints,
    )
    return _conformation_types(st, conformations, conformation_files)
end
