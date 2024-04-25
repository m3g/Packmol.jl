@kwdef mutable struct PackmolSystem{N,T}
    tolerance::T = 2.0
    structure_types::Vector{StructureType{N,T}} = StructureType{N,T}[]
    input_file::String = "none"
    output_file::String = "packmol_output.pdb"
    tolerance_precision::T = 1e-2
    constraint_precision::T = 1e-2
    max_iter::Int = 1000
end