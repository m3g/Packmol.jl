#=

    StructureType

Structure that contains the input data for a structure block in the input file. 

=#
@kwdef struct StructureType{D,T}
    filename::String
    natoms::Int
    atoms::Vector{Atom{Nothing}}
    number_of_molecules::Int
    fixed::FixedMoleculeData{D,T} = zero(FixedMoleculeData{D,T})
    reference_coordinates::Vector{SVector{D,T}}
    radii::Vector{T} = T[]
    # -1 (unset, Fortran Packmol's own sentinel) is resolved at output time
    # (write_output.jl) to 0 for a single-residue template or 1 for a
    # multi-residue one — matching Fortran's `setrnum` auto-detection, since
    # a literal default of 1 (verbatim template residue number) would give
    # every copy of a single-residue solvent molecule the same residue
    # number.
    residue_numbering::Int = -1
    # `nothing` (Fortran's `chain` unset, "#") means auto-assign a chain letter
    # per molecule at output time (see write_output.jl); an explicit `Char`
    # (from the `chain <c>` keyword) fixes every molecule of this type to
    # that one chain identifier instead.
    chain::Union{Nothing,Char} = nothing
    changechains::Bool = false
    # Segment identifier (`segid` keyword) written to every atom of this type
    # in the output (PDB columns 73-76); `nothing` keeps the template's own.
    segid::Union{Nothing,String} = nothing
    # `constrain_rotation`: bounds on the angle between an axis and a vector of the
    # molecule (see constraints/rotations.jl)
    rotation_constraints::Vector{RotationConstraint{D,T}} = RotationConstraint{D,T}[]
    # Per-structure-type restart (Fortran's `restart_from`/`restart_to` given
    # inside a `structure ... end structure` block): `nothing` means unset.
    # See `PackmolSystem.restart_from`/`restart_to` for the whole-system
    # version, and its docstring for the accepted `restart_from` sources.
    restart_from::Union{Nothing,String,Vector{Atom{Nothing}}} = nothing
    restart_to::Union{Nothing,String} = nothing
    # Bonds of the template (from the CONECT records of the structure file, read by
    # PDBToolsBase): `connectivity[i]` are the indices of the atoms bonded to atom
    # `i`. `nothing` if not available, or disabled with `connect no`.
    connectivity::Union{Nothing,Vector{Vector{Int32}}} = nothing
    # AnyConstraint{T} (a closed Union of every concrete constraint type, not
    # the abstract Constraint) so constraint dispatch in the hot packing loop
    # is allocation-free — see constraints/constraint_types.jl.
    constraints::Vector{AnyConstraint{T}}
    atom_constraints::Vector{Vector{Int}} = Vector{Int}[]
    move_bad_molecules::Symbol = :low_density_region
    # Conformations: a molecule given in several conformations (a structure file
    # with several models, or a vector of files) is packed as consecutive structure
    # types, one per conformation; `conformation = (i, n)` marks this one as the
    # i-th of a set of n (`(1, 1)` for a molecule with a single conformation).
    conformation::Tuple{Int,Int} = (1, 1)
end

function _show(s::StructureType{D,T}) where {D,T}
    chomp(
        """
        StructureType{$D,$T}
            filename: $(basename(s.filename))
            natoms: $(s.natoms)
            number_of_molecules: $(s.number_of_molecules)
            fixed: $(s.fixed.fixed)
            number of constraints: $(length(s.constraints))$(s.conformation[2] > 1 ? "\n    conformation: $(s.conformation[1]) of $(s.conformation[2])" : "")
        """
    )
end
Base.show(io::IO, ::MIME"text/plain", s::StructureType) = print(io, _show(s))

# One-line summary of a structure type (in the `show` of a `PackmolSystem`), or of a
# set of conformations of a molecule (`sts`, consecutive structure types), as a single
# structure type
_summary(s::StructureType) = _summary([s])
function _summary(sts::AbstractVector{<:StructureType})
    s = first(sts)
    n = length(sts)
    return string(
        "\"", basename(s.filename), "\"", n > 1 ? " - $n conformations" : "",
        " - n_atoms = ", s.natoms, " - n_mols = ", sum(st.number_of_molecules for st in sts),
        s.fixed.fixed ? " - fixed" : "", " - n_constraints = ", length(s.constraints),
    )
end

# The structure types grouped by molecule: each set of conformations of a molecule is
# one group (ranges of indices of `structure_types`)
function _molecule_groups(structure_types::AbstractVector{<:StructureType})
    groups = UnitRange{Int}[]
    for (i, st) in enumerate(structure_types)
        if st.conformation[1] == 1 || isempty(groups)
            push!(groups, i:i)
        else
            groups[end] = first(groups[end]):i
        end
    end
    return groups
end

function Base.show(io::IO, ::MIME"text/plain", v::AbstractVector{<:StructureType})
    print(io, chomp("""
        Vector{StructureType} with $(length(v)) structure(s).
    """))
end

#
# Apply the `fixed`/`center` transformation to `reference_coordinates` in place:
# if `fixed.fixed`, center the molecule (per `center`: `:mass`, `:geometric`, or
# neither) and then move/rotate it to `fixed.position`; otherwise, `center`
# must be `:none`. Shared by `read_structure_data` and `structure_type`.
#
function _apply_fixed_center!(
    reference_coordinates::Vector{SVector{D,T}}, atoms, fixed::FixedMoleculeData{D,T}, center::Symbol,
) where {D,T}
    if fixed.fixed
        if center == :mass
            # centerofmass: use mass-weighted center of mass from PDBToolsBase
            com = center_of_mass(atoms)
            cm_vec = SVector{D,T}(ntuple(i -> com[i], D))
            reference_coordinates .-= Ref(cm_vec)
            R = eulermat(fixed.position.angles)
            for i in eachindex(reference_coordinates)
                reference_coordinates[i] = R * reference_coordinates[i]
            end
            reference_coordinates .+= Ref(fixed.position.cm)
        elseif center == :geometric
            move!(reference_coordinates, fixed.position)
        else
            # Without `center`/`centerofmass`, the molecule is rotated about the origin of
            # its coordinates and then translated, as in Fortran Packmol
            R = eulermat(fixed.position.angles)
            for i in eachindex(reference_coordinates)
                reference_coordinates[i] = R * reference_coordinates[i] + fixed.position.cm
            end
        end
    elseif center != :none
        throw(ArgumentError("option `center`/`centerofmass` cannot be set without a fixed position"))
    end
    return reference_coordinates
end

# The position of an atom, as an `SVector{D,T}`
_atom_position(atom::Atom, ::Val{D}, ::Type{T}) where {D,T} = SVector{D,T}(ntuple(i -> (atom.x, atom.y, atom.z)[i], D))

"""
    structure_type(filename; number, constraints=Constraint[], tolerance=2.0, radius=nothing,
                    fixed=nothing, center=:none, residue_numbering=-1, D=3, T=Float64)

Build a `StructureType` for `number` copies of the molecule read from the PDB file
`filename`, without going through an input file. Every constraint in `constraints`
applies to the whole molecule.

`atom_constraints` gives constraints that apply only to some atoms of the molecule (as
the `atoms ... end atoms` blocks of an input file): a vector of `indices => constraints`
pairs, with the indices of the atoms in the structure file and one constraint or a
vector of them. For example, `atom_constraints = [[31, 32] => BelowPlane([0,0,1], 2.0)]`.

`filename` can also be a vector of structure files, or a file with several models, for
several conformations of the molecule (they must have the same atoms, in the same
order): a vector of structure types is then returned, one per conformation, among
which `number` is divided as evenly as possible, and which `PackmolSystem` takes as a
single structure type (e.g. `PackmolSystem([lipids, water]; ...)`).

`radius` sets a uniform atom radius; if not given, it defaults to `tolerance/2` (so
that two touching atoms' radii sum to `tolerance`, matching the input file default).

`fixed`, if given, is a `(cm, angles)` pair of 3-vectors that fixes the molecule at
that position/orientation instead of packing it; `center` (`:none`, `:geometric`, or
`:mass`) then controls how the molecule is centered before that transformation is
applied — equivalent to the input file's `fixed`/`center`/`centerofmass` keywords.

# Example

```julia-repl
julia> using Packmol

julia> st = structure_type("water.pdb"; number=1000, constraints=[InsideBox([0,0,0],[40,40,40])])
```
"""
function structure_type(filename::Union{String,AbstractVector{String}};
    number::Int,
    constraints::Vector{<:Constraint} = Constraint[],
    atom_constraints::AbstractVector{<:Pair} = Pair[],
    tolerance::Real = 2.0,
    radius::Union{Nothing,Real} = nothing,
    fixed::Union{Nothing,Tuple} = nothing,
    center::Symbol = :none,
    residue_numbering::Int = -1,
    chain::Union{Nothing,Char} = nothing,
    changechains::Bool = false,
    segid::Union{Nothing,String} = nothing,
    # e.g. `constrain_rotation=Dict(:z => (0, 10))` keeps the angle between the z axis and
    # the longest axis of the molecule within 0° ± 10°, and `Dict(:z => (0, 10, [1, 2], [30]))`
    # the angle between the z axis and the vector from the center of atoms 1 and 2 to atom 30
    constrain_rotation = Dict{Symbol,Tuple{Float64,Float64}}(),
    restart_from::Union{Nothing,String,AbstractVector{<:Atom}} = nothing,
    restart_to::Union{Nothing,String} = nothing,
    connect::Bool = true,
    filetype::String = "pdb",
    D::Int = 3,
    T::DataType = Float64,
)
    changechains && !isnothing(chain) &&
        error("'changechains' and 'chain' are not compatible for a single structure.")
    conformations, conformation_files, connectivity = _read_conformations(filename, filetype)
    atoms = conformations[1]
    filename = conformation_files[1]
    connect || (connectivity = nothing)
    natoms = length(atoms)
    reference_coordinates = [_atom_position(at, Val(D), T) for at in atoms]
    radii = fill(T(something(radius, tolerance / 2)), natoms)
    per_atom_constraints = atom_constraints
    constraints = Constraint[constraints...]
    atom_constraints = [collect(1:length(constraints)) for _ in 1:natoms]
    for (indices, cs) in per_atom_constraints
        for i in indices
            (1 <= i <= natoms) || throw(ArgumentError("""\n
                atom_constraints: atom index $i out of range (the molecule has $natoms atoms).

            """))
        end
        for c in (cs isa Constraint ? (cs,) : cs)
            push!(constraints, c)
            for i in indices
                push!(atom_constraints[i], length(constraints))
            end
        end
    end

    fixed_data = if isnothing(fixed)
        zero(FixedMoleculeData{D,T})
    else
        cm = SVector{D,T}(fixed[1]...)
        angles = SVector{D,T}(fixed[2]...)
        FixedMoleculeData(true, MoleculePosition(cm, angles))
    end
    _apply_fixed_center!(reference_coordinates, atoms, fixed_data, center)

    rotation_constraints = RotationConstraint{D,T}[]
    for (axis, spec) in constrain_rotation
        length(spec) in (2, 4) || throw(ArgumentError("""\n
            constrain_rotation: give (center, halfwidth) or (center, halfwidth, from_atoms, to_atoms) for each axis.

        """))
        from, to = length(spec) == 4 ? (spec[3], spec[4]) : (Int[], Int[])
        push!(rotation_constraints, RotationConstraint{D,T}(axis, spec[1], spec[2]; from, to))
    end
    _check_rotation_constraints(rotation_constraints, fixed_data, natoms)

    st = StructureType{D,T}(;
        filename, natoms, atoms, number_of_molecules=number, rotation_constraints, restart_from=_restart_source(restart_from), restart_to,
        fixed=fixed_data, reference_coordinates, radii,
        residue_numbering, chain, changechains, segid=_check_segid(segid), connectivity, constraints, atom_constraints,
    )
    types = _conformation_types(st, conformations, conformation_files)
    return length(types) == 1 ? only(types) : types
end

#
# The conformations of a molecule (see `_conformations`), from one structure file or
# from several, with the file of each conformation, and the connectivity of the first.
#
function _read_conformations(filenames::Union{String,AbstractVector{String}}, filetype::String)
    filenames isa String && (filenames = [filenames])
    isempty(filenames) && throw(ArgumentError("""\n
        No structure file given.

    """))
    conformations = Vector{Atom{Nothing}}[]
    conformation_files = String[]
    connectivity = nothing
    for (i, file) in enumerate(filenames)
        atoms, conn = _read_structure_file(file, filetype)
        i == 1 && (connectivity = conn)
        for conf in _conformations(atoms)
            push!(conformations, conf)
            push!(conformation_files, file)
        end
    end
    # The bonds of the first conformation (the CONECT records of a file with several
    # models refer to the atoms of the first one)
    n = length(conformations[1])
    if !isnothing(connectivity) && length(connectivity) != n
        connectivity = length(connectivity) > n && all(all(<=(n), c) for c in @view(connectivity[1:n])) ?
            connectivity[1:n] : nothing
    end
    return conformations, conformation_files, connectivity
end

read_structure_data(input_file_block::AbstractString, args...; kargs...) = 
    read_structure_data(IOBuffer(input_file_block), args...; kargs...)

read_structure_data(input_file_block::IOBuffer, tolerance; D::Int=3, T::DataType=Float64, path::String="", filetype::String="pdb") =
    read_structure_data(input_file_block, tolerance, Val(D), T; path, filetype)

#
# Atoms and connectivity of a structure file. The format follows the extension
# of the file (`.pdb`, `.xyz`, `.cif`/`.mmcif`) or, for other extensions, the
# `filetype` of the input (`pdb` or `xyz`, as in Fortran Packmol).
#
function _read_structure_file(filename::String, filetype::String)
    structure = if _file_format(filename, filetype) == "xyz"
        read_xyz(filename)
    elseif _file_format(filename, filetype) == "mmcif"
        read_mmcif(filename)
    else
        read_pdb(filename)
    end
    return Vector{Atom{Nothing}}(structure.atoms), structure.connectivity
end

#
# The conformations of the molecule of a structure file: the atoms of each of its
# models (`MODEL` records of a PDB file), or all its atoms if it has a single model.
#
function _conformations(atoms::Vector{Atom{Nothing}})
    models = Int32[]
    for at in atoms
        at.model in models || push!(models, at.model)
    end
    length(models) == 1 && return [atoms]
    conformations = [Atom{Nothing}[] for _ in models]
    for at in atoms
        at = copy(at)
        k = findfirst(==(at.model), models)
        at.model = 1
        push!(conformations[k], at)
    end
    return conformations
end

#
# A molecule given in several conformations is packed as consecutive structure
# types, one per conformation, marked by their `conformation` field: copies of `st`
# (built from the first conformation) with the coordinates of each conformation,
# and the number of molecules split as evenly as possible among them (conformations
# left without molecules are dropped). Each copy has its own constraint vectors,
# since constraints may be added to each structure type later (e.g. the periodic
# cell confinement).
#
function _conformation_types(
    st::StructureType{D,T}, conformations::Vector{Vector{Atom{Nothing}}}, filenames::Vector{String},
) where {D,T}
    nconf = length(conformations)
    nconf == 1 && return [st]
    st.fixed.fixed && throw(ArgumentError("""\n
        Structure $(basename(st.filename)): a fixed molecule cannot have several conformations (models).

    """))
    (isnothing(st.restart_from) && isnothing(st.restart_to)) || throw(ArgumentError("""\n
        Structure $(basename(st.filename)): restart_from and restart_to are not supported for a molecule
        with several conformations (models) in the structure block; use them for the whole system.

    """))
    for (k, conf) in enumerate(conformations)
        if length(conf) != st.natoms || any(conf[i].name != st.atoms[i].name for i in 1:st.natoms)
            throw(ArgumentError("""\n
                Structure $(basename(st.filename)): conformation $k ($(basename(filenames[k]))) does not have the
                same atoms, in the same order, as the first one.

            """))
        end
    end
    q, r = divrem(st.number_of_molecules, nconf)
    keep = [k for k in 1:nconf if q + (k <= r) > 0]
    types = StructureType{D,T}[]
    for (i, k) in enumerate(keep)
        atoms = conformations[k]
        push!(types, StructureType{D,T}(;
            filename=filenames[k], natoms=st.natoms, atoms,
            number_of_molecules=q + (k <= r), fixed=st.fixed,
            reference_coordinates=[_atom_position(at, Val(D), T) for at in atoms],
            radii=copy(st.radii), residue_numbering=st.residue_numbering, chain=st.chain,
            changechains=st.changechains, segid=st.segid, rotation_constraints=copy(st.rotation_constraints),
            restart_from=nothing, restart_to=nothing, connectivity=st.connectivity,
            constraints=copy(st.constraints), atom_constraints=[copy(ac) for ac in st.atom_constraints],
            move_bad_molecules=st.move_bad_molecules, conformation=(i, length(keep)),
        ))
    end
    return types
end

function _file_format(filename::AbstractString, filetype::String)
    ext = lowercase(splitext(filename)[2])
    ext == ".pdb" && return "pdb"
    ext == ".xyz" && return "xyz"
    (ext == ".cif" || ext == ".mmcif") && return "mmcif"
    (filetype == "pdb" || filetype == "xyz") && return filetype
    throw(ArgumentError("""\n
        filetype $filetype is not supported: use pdb or xyz.

    """))
end

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
    input_file_block, rotation_constraints = _read_rotation_constraints(input_file_block, Val(D), T)
    restart_from = nothing
    restart_to = nothing
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
        elseif keyword == "restart_from"
            restart_from = joinpath(path, expanduser(String(values[1])))
        elseif keyword == "restart_to"
            restart_to = joinpath(path, expanduser(String(values[1])))
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
        residue_numbering, chain, changechains, segid, rotation_constraints, restart_from, restart_to,
        connectivity=connect ? connectivity : nothing,
        constraints, atom_constraints,
    )
    return _conformation_types(st, conformations, conformation_files)
end
