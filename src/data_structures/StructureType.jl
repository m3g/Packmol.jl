#=

    StructureType

Structure that contains the input data for a structure block in the input file. 

=#
@kwdef struct StructureType{D,T}
    filename::String
    natoms::Int
    atoms::Vector{Atom}
    number_of_molecules::Int
    fixed::Bool = false
    reference_coordinates::Vector{SVector{D,T}}
    radii::Vector{T} = T[]
    residue_numbering::Int = 1
    connect::Vector{Vector{Int}} = Vector{Int}[]
    constraints::Vector{<:Constraint}
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
            fixed: $(s.fixed)
            number of constraints: $(length(s.constraints))
        """
    )
end
Base.show(io::IO, ::MIME"text/plain", s::StructureType) = print(io, _show(s))

#=
    read_structure_data(input_file_block::IOBuffer, tolerance; T=Float64, D=3)
    read_structure_data(input_file_block::AbstractString, args... kargs..)

Reads the structure data from a structure block of the input file. Requires the `tolerance`
parameter to set atom radii by default. The function returns a `StructureType` object. 

The type `T` defines the number type (Float32, Float64, etc.) and `D` the number of dimensions (2 or 3).

=#
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
        :constraints => Constraint[]
    )
    # Read basic structure data first
    seekstart(input_file_block)
    atoms_block = false
    iconstraint = 0
    for line in eachline(input_file_block)
        keyword, values... = split(line)
        if keyword == "atoms"
            atoms_block = true
        end
        if keyword == "end" && atoms_block
            atoms_block = false
        end
        if keyword == "structure"
            filename = joinpath(path, string(values[1]))
            atoms = readPDB(filename)
            structure_data[:filename] = filename
            structure_data[:natoms] = length(atoms)
            structure_data[:atoms] = atoms
            structure_data[:reference_coordinates] = coor(atoms)
            structure_data[:radii] = fill(T(tolerance), structure_data[:natoms])
            structure_data[:atom_constraints] = [Int[] for _ in 1:structure_data[:natoms]]
        elseif keyword == "number"
            structure_data[:number_of_molecules] = parse(Int, values[1]) 
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

@testitem "read_structure_data" begin
    file = Packmol.src_dir*"/../test/data/water.pdb"
    tolerance = 2.0

    input_file_block = """
    structure $file        
        number 1000
        inside box 0. 0. 0. 40. 40. 40.
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance)
    @test isfile(s.filename)
    @test length(s.atoms) == 3
    @test s.number_of_molecules == 1000
    @test length(s.reference_coordinates) == 3
    @test length(s.constraints) == 1
    @test s.constraints[1] == Box{Inside, Float64}([0.0, 0.0, 0.0], [40.0, 40.0, 40.0], 5.0)
    @test s.fixed == false 
    @test s.atom_constraints == [[1], [1], [1]]
    @test s.radii == [2.0, 2.0, 2.0]

    input_file_block = """
    structure $file        
        number 1000
        inside box 0. 0. 0. 40. 40. 40.
        outside sphere 0. 0. 0. 10.
        radius 1.0
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance)
    @test s.radii == [1.0, 1.0, 1.0]
    @test s.atom_constraints == [[1, 2], [1, 2], [1, 2]]
    @test s.constraints[1] == Box{Inside, Float64}([0.0, 0.0, 0.0], [40.0, 40.0, 40.0], 5.0)
    @test s.constraints[2] == Sphere{Outside, Float64}([0.0, 0.0, 0.0], 10.0, 5.0)

    input_file_block = """
    structure $file        
        number 1000
        inside box 0. 0. 0. 40. 40. 40.
        atoms 1 3
            outside sphere 0. 0. 0. 10.
        end
        atoms 1 2
            radius 1.0
        end
    end structure
    """
    s = Packmol.read_structure_data(input_file_block, tolerance)
    @test s.radii == [1.0, 1.0, 2.0]
    @test s.atom_constraints == [[1, 2], [1], [1, 2]]
    @test s.constraints[1] == Box{Inside, Float64}([0.0, 0.0, 0.0], [40.0, 40.0, 40.0], 5.0)
    @test s.constraints[2] == Sphere{Outside, Float64}([0.0, 0.0, 0.0], 10.0, 5.0)

    s = Packmol.read_structure_data(input_file_block, tolerance; T = Float32)
    @test s.radii == Float32[1.0, 1.0, 2.0]
    @test s.constraints[1] == Box{Inside, Float32}([0.0, 0.0, 0.0], [40.0, 40.0, 40.0], 5.0)
    @test s.constraints[2] == Sphere{Outside, Float32}([0.0, 0.0, 0.0], 10.0, 5.0)
end


