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
    # Per-axis (index 1=x, 2=y, 3=z — matching `eulermat`'s angle order, see
    # rigid_body.jl) bound on that molecule's rotation angle, as
    # (lower, upper) in radians; `nothing` means unconstrained. Set by
    # `constrain_rotation <axis> <center_deg> <halfwidth_deg>`, one line per
    # constrained axis. Applied as hard bounds on the optimization variables
    # in the main packing loop (matching Fortran Packmol's own `pgencan`,
    # which sets GENCAN's `l`/`u` this way) and as the sampling range for
    # that axis during initial random placement.
    rotation_bounds::Vector{Union{Nothing,Tuple{T,T}}} = fill(nothing, D)
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
end

function _show(s::StructureType{D,T}) where {D,T}
    chomp(
        """
        StructureType{$D,$T}
            filename: $(basename(s.filename))
            natoms: $(s.natoms)
            number_of_molecules: $(s.number_of_molecules)
            fixed: $(s.fixed.fixed)
            number of constraints: $(length(s.constraints))
        """
    )
end
Base.show(io::IO, ::MIME"text/plain", s::StructureType) = print(io, _show(s))

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
            cm = mean(reference_coordinates)
            move!(reference_coordinates, fixed.position)
            reference_coordinates .+= Ref(cm)
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
applies to the whole molecule (for per-atom constraints, build a `StructureType`
directly with its keyword constructor).

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
function structure_type(filename::String;
    number::Int,
    constraints::Vector{<:Constraint} = Constraint[],
    tolerance::Real = 2.0,
    radius::Union{Nothing,Real} = nothing,
    fixed::Union{Nothing,Tuple} = nothing,
    center::Symbol = :none,
    residue_numbering::Int = -1,
    chain::Union{Nothing,Char} = nothing,
    changechains::Bool = false,
    segid::Union{Nothing,String} = nothing,
    # e.g. `constrain_rotation=Dict(:x => (45, 10))` bounds the x-axis
    # rotation angle to 45° ± 10°; see `StructureType.rotation_bounds`.
    constrain_rotation::AbstractDict{Symbol} = Dict{Symbol,Tuple{Float64,Float64}}(),
    restart_from::Union{Nothing,String,AbstractVector{<:Atom}} = nothing,
    restart_to::Union{Nothing,String} = nothing,
    connect::Bool = true,
    filetype::String = "pdb",
    D::Int = 3,
    T::DataType = Float64,
)
    changechains && !isnothing(chain) &&
        error("'changechains' and 'chain' are not compatible for a single structure.")
    atoms, connectivity = _read_structure_file(filename, filetype)
    connect || (connectivity = nothing)
    natoms = length(atoms)
    reference_coordinates = [_atom_position(at, Val(D), T) for at in atoms]
    radii = fill(T(something(radius, tolerance / 2)), natoms)
    atom_constraints = [collect(1:length(constraints)) for _ in 1:natoms]

    fixed_data = if isnothing(fixed)
        zero(FixedMoleculeData{D,T})
    else
        cm = SVector{D,T}(fixed[1]...)
        angles = SVector{D,T}(fixed[2]...)
        FixedMoleculeData(true, MoleculePosition(cm, angles))
    end
    _apply_fixed_center!(reference_coordinates, atoms, fixed_data, center)

    rotation_bounds = Union{Nothing,Tuple{T,T}}[nothing for _ in 1:D]
    axis_index = Dict(:x => 1, :y => 2, :z => 3)
    for (axis, (center_deg, halfwidth_deg)) in constrain_rotation
        haskey(axis_index, axis) || error("constrain_rotation axis must be :x, :y, or :z, got $(repr(axis))")
        idx = axis_index[axis]
        idx > D && error("constrain_rotation axis $axis is out of range for D=$D")
        hw = abs(T(halfwidth_deg)) * T(π) / 180
        c = T(center_deg) * T(π) / 180
        rotation_bounds[idx] = (c - hw, c + hw)
    end

    return StructureType{D,T}(;
        filename, natoms, atoms, number_of_molecules=number, rotation_bounds, restart_from=_restart_source(restart_from), restart_to,
        fixed=fixed_data, reference_coordinates, radii,
        residue_numbering, chain, changechains, segid=_check_segid(segid), connectivity, constraints, atom_constraints,
    )
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

function read_structure_data(input_file_block::IOBuffer, tolerance, ::Val{D}, ::Type{T}; path::String="", filetype::String="pdb") where {D,T}
    filename = ""
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
    rotation_bounds = Union{Nothing,Tuple{T,T}}[nothing for _ in 1:D]
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
            atoms, connectivity = _read_structure_file(filename, filetype)
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
        elseif keyword == "constrain_rotation"
            axis, center_str, halfwidth_str = values[1], values[2], values[3]
            idx = axis == "x" ? 1 : axis == "y" ? 2 : axis == "z" ? 3 :
                error("constrain_rotation axis must be x, y, or z, got '$axis'")
            idx > D && error("constrain_rotation axis '$axis' is out of range for D=$D")
            hw = abs(parse(T, halfwidth_str)) * T(π) / 180
            c = parse(T, center_str) * T(π) / 180
            rotation_bounds[idx] = (c - hw, c + hw)
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
    return StructureType{D,T}(;
        filename, natoms=length(atoms), atoms, number_of_molecules, fixed, reference_coordinates, radii,
        residue_numbering, chain, changechains, segid, rotation_bounds, restart_from, restart_to,
        connectivity=connect ? connectivity : nothing,
        constraints, atom_constraints,
    )
end
