#=

    StructureType

Structure that contains the input data for a structure block in the input file. 

=#
@kwdef struct StructureType{D,T}
    filename::String
    natoms::Int
    atoms::Vector{Atom}
    number_of_molecules::Int
    fixed::FixedMoleculeData{D,T} = zero(FixedMoleculeData{D,T})
    reference_coordinates::Vector{SVector{D,T}}
    radii::Vector{T} = T[]
    residue_numbering::Int = 1
    connect::Vector{Vector{Int}} = Vector{Int}[]
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
            # centerofmass: use mass-weighted center of mass from PDBTools
            cm_vec = SVector{D,T}(center_of_mass(atoms)[1:D]...)
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

"""
    structure_type(filename; number, constraints=Constraint[], tolerance=2.0, radius=nothing,
                    fixed=nothing, center=:none, residue_numbering=1, D=3, T=Float64)

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
    residue_numbering::Int = 1,
    D::Int = 3,
    T::DataType = Float64,
)
    atoms = read_pdb(filename)
    natoms = length(atoms)
    reference_coordinates = [SVector{D,T}(c[1:D]...) for c in coor(atoms)]
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

    return StructureType{D,T}(;
        filename, natoms, atoms, number_of_molecules=number,
        fixed=fixed_data, reference_coordinates, radii,
        residue_numbering, constraints, atom_constraints,
    )
end

read_structure_data(input_file_block::AbstractString, args...; kargs...) = 
    read_structure_data(IOBuffer(input_file_block), args...; kargs...)

function read_structure_data(input_file_block::IOBuffer, tolerance; 
    D::Int=3, T::DataType=Float64, path::String="",
)
    constraint_placements = first.(split.(keys(parse_constraint)))
    structure_data = Dict{Symbol,Any}(
        :filename => nothing,
        :natoms => nothing,
        :atoms => nothing,
        :number_of_molecules => nothing,
        :reference_coordinates => nothing,
        :constraints => Constraint[],
        :fixed => zero(FixedMoleculeData{D,T}),
        :center => :none,  # :none, :geometric, or :mass
    )
    # Read basic structure data first
    seekstart(input_file_block)
    atoms_block = false
    iconstraint = 0
    seekstart(input_file_block)
    for line in eachline(input_file_block)
        keyword, values... = split(line)
        if keyword == "atoms"
            atoms_block = true
        end
        if atoms_block && keyword == "end"
            atoms_block = false
        end
        if keyword == "structure"
            filename = joinpath(path, string(values[1]))
            atoms = read_pdb(filename)
            structure_data[:filename] = filename
            structure_data[:natoms] = length(atoms)
            structure_data[:atoms] = atoms
            structure_data[:reference_coordinates] = [SVector{D,T}(c[1:D]...) for c in coor(atoms)]
            structure_data[:radii] = fill(T(tolerance/2), structure_data[:natoms])
            structure_data[:atom_constraints] = [Int[] for _ in 1:structure_data[:natoms]]
        elseif keyword == "number"
            structure_data[:number_of_molecules] = parse(Int, values[1]) 
        elseif keyword == "fixed"
            vals = parse.(T, values)
            cm = SVector{D,T}(vals[1:D]...)
            angles = SVector{D,T}(vals[D+1:2*D]...)
            structure_data[:fixed] = FixedMoleculeData(true, MoleculePosition(cm, angles))
        elseif keyword == "center"
            structure_data[:center] = :geometric
        elseif keyword == "centerofmass"
            structure_data[:center] = :mass
        elseif keyword == "resnumbers"
            structure_data[:residue_numbering] = parse(Int, values[1])
        elseif keyword in constraint_placements
            iconstraint += 1
            push!(structure_data[:constraints], parse_constraint["$keyword $(values[1])"](structure_data, values[2:end]; T=T))
            if !atoms_block
                for iatom in 1:structure_data[:natoms]
                    push!(structure_data[:atom_constraints][iatom], iconstraint)
                end
            end
        elseif !atoms_block && keyword == "radius"
            structure_data[:radii] .= parse(T, values[1])
        end
    end
    # If molecule is fixed, apply transformation to obtain the reference coordinates
    _apply_fixed_center!(
        structure_data[:reference_coordinates], structure_data[:atoms],
        structure_data[:fixed], structure_data[:center],
    )
    pop!(structure_data, :center)
    #
    # Read custom atom radii and set constraints, according to atom blocks
    #
    atoms_block = false
    iconstraint = 0
    seekstart(input_file_block)
    local atoms_in_block
    for line in eachline(input_file_block)
        keyword, values... = split(line)
        if keyword == "atoms" 
            atoms_block = true
            atoms_in_block = parse.(Int, values)
        end
        if keyword == "end" && atoms_block
            atoms_block = false
        end
        if atoms_block && keyword == "radius"
            radius = parse(T, values[1])
            structure_data[:radii][atoms_in_block] .= radius
        end
        if keyword in constraint_placements
            iconstraint += 1
            if atoms_block 
                for iatom in atoms_in_block
                    push!(structure_data[:atom_constraints][iatom], iconstraint)
                end
            end
        end
    end
    return StructureType{D,T}(;structure_data...)
end
