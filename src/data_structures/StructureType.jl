@kwdef struct StructureType{N,T}
    filename::String
    number_of_molecules::Int
    fixed::Bool = false
    reference_coordinates::Vector{SVector{N,T}}
    radii::Vector{T} = [ zero(T) for _ in eachindex(reference_coordinates) ]
    residue_numbering::Int = 1
    connect::Vector{Vector{Int}} = Vector{Int}[]
    constraints::Vector{<:Constraint}
    move_bad_molecules::Symbol = :low_density_region
end

function read_structure_data(input_file_block::String; T=Float64, N=3)
    constraint_placements = first.(split.(keys(parse_constraint)))
    structure_data = Dict{Symbol,Any}(
        :filename => nothing,
        :number_of_molecules => nothing,
        :reference_coordinates => nothing,
        :constraints => Constraint[]
    )
    for line in eachline(IOBuffer(input_file_block))
        data = split(line)
        if data[1] == "structure"
            filename = String(data[2])
            structure_data[:filename] = filename;
            structure_data[:reference_coordinates] = coor(readPDB(filename))
        elseif data[1] == "number"
            structure_data[:number_of_molecules] = parse(Int, data[2]) 
        elseif data[1] in constraint_placements
            push!(structure_data[:constraints], parse_constraint["$(data[1]) $(data[2])"](structure_data, data; T=T, N=N))
        elseif data[1] == "end"
            break
        end
    end
    return StructureType{N,T}(;structure_data...)
end

@testitem "read_structure_data" begin
    input_file_block = """
    structure water.pdb
        number 1000
        inside box 0. 0. 0. 40. 40. 40.
    end structure
    """

end


