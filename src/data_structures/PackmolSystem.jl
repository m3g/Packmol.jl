@kwdef mutable struct PackmolSystem{N,T}
    input_file::String
    output_file::String = "packmol_output.pdb"
    tolerance::T = 2.0
    structure_types::Vector{StructureType{N,T}} = StructureType{N,T}[]
    tolerance_precision::T = 1e-2
    constraint_precision::T = 1e-2
    max_iter::Int = 1000
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
end


#=

    read_packmol_input

Reads the packmol input file to create a `PackmolSystem` object.

=#
function read_packmol_input(input_file::String; T=Float64, N=3)
    packmol_input_data = Dict{Symbol,Any}(
        :input_file => input_file,
    )
    for line in eachline(input_file)
        keyword, values... = split(line)
        if keyword == "output"
            packmol_input_data[:output_file] = values[1]
        elseif keyword == "tolerance"
            packmol_input_data[:tolerance] = parse(T, values[1])
        elseif keyword == "tolerance_precision"
            packmol_input_data[:tolerance_precision] = parse(T, values[1])
        elseif keyword == "constraint_precision"
            packmol_input_data[:tolerance_precision] = parse(T, values[1])
        elseif keyword == "precision" 
            packmol_input_data[:tolerance_precision] = parse(T, values[1])
            packmol_input_data[:constraint_precision] = parse(T, values[1])
        elseif keyword == "max_iter"
            packmol_input_data[:max_iter] = parse(Int, values[1])
        elseif keyword == "add_amber_ter"
            packmol_input_data[:add_amber_ter] = true
        elseif keyword == "amber_ter_preserve"
            packmol_input_data[:amber_ter_preserve] = true
        elseif keyword == "add_box_sides"
            packmol_input_data[:add_amber_ter] = true
        elseif keyword == "connect"
            packmol_input_data[:connect] = values[1] == "yes" ? true : false
        elseif keyword == "randominitialpoint"
            packmol_input_data[:random_initial_point] = true
        elseif keyword == "seed"
            packmol_input_data[:seed] = parse(Int, values[1])
        elseif keyword == "avoid_overlap"
            packmol_input_data[:avoid_overlap] = values[1] == "yes" ? true : false
        elseif keyword == "writeout"
            packmol_input_data[:writeout] = parse(Int,values[1])
        elseif keyword == "writebad"
            packmol_input_data[:writebad] = true
        elseif keyword == "optim_print_level"
            packmol_input_data[:optim_print_level] = parse(Int,values[1])
        elseif keyword == "chkgrad"
            packmol_input_data[:chkgrad] = true
        end
    end
    return PackmolSystem{N,T}(;packmol_input_data...)
end