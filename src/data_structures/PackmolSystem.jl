struct AtomData{D,T}
    molecule_index::Int
    structure_type_index::Int
    radius::T
    constraints::Vector{Constraint}
end

# MoleculePosition: This is a central data structure that 
# contains the center of mass and the rotation angles for each molecule.
# This data structure is used to build the array that contains 
# rotation angles and center of mass for each molecule, which are 
# the variables of the optimization problem. The data structure can,
# and will, be reinterpreted as a linear vector, to conform with the 
# interface of the optimization rotations, using:
#
#     reinterpret(Float64, molecule_positions)
#  
# The resulting vector contains, in order, the center of mass and the
# rotation angles for each molecule. The same data structure will
# be used to store the gradient of the objective function relative to
# the rotations and translation of the rigid-body molecules.
struct MoleculePosition{N,T}
    cm::SVector{N,T}
    angles::SVector{N,T}
end
Base.copy(x::MoleculePosition) = MoleculePosition(x.cm, x.angles)
import Base: + 
+(x::MoleculePosition, y::MoleculePosition) = MoleculePosition(x.cm + y.cm, x.angles + y.angles)
Base.zero(::Type{MoleculePosition{N,T}}) where {N,T} = MoleculePosition(zero(SVector{N,T}), zero(SVector{N,T}))

@kwdef mutable struct PackmolSystem{D,T}
    filetype::String
    input_file::String
    output_file::String
    tolerance::T = 2.0
    structure_types::Vector{StructureType{D,T}} = StructureType{D,T}[]
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
    atoms::Vector{AtomData{D,T}} = AtomData{D,T}[]
    nmols::Int = 0
    molecule_positions::MoleculePosition{D,T} = MoleculePosition{D,T}[]
end

function _indent(s::AbstractString; n=4)
    idented_str = IOBuffer()
    for line in eachline(IOBuffer(s))
        println(idented_str, repeat(" ", n) * line)
    end
    return String(take!(idented_str))
end

function Base.show(io::IO, mime::MIME"text/plain", sys::PackmolSystem{D,T}) where {D,T}
    printstyled(io, "PackmolSystem{$D,$T}"; bold=true, color=:blue)
    if length(sys.input_file) > 0 
        printstyled(" - read from: $(basename(sys.input_file))"; color=:blue) 
    end
    println(io)
    printstyled(io, _indent("Structure types:\n"); bold=true)
    for (i,st) in enumerate(sys.structure_types)
        print(io, _indent("$i."*_show(st)))
    end
    printstyled(io, _indent("Options:\n"); bold=true)
    for field in fieldnames(PackmolSystem)
        if !(field in (:structure_types, :input_file))
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

_check_movefrac(x) = 0.0 <= x <= 1.0 ? x : throw(ArgumentError("movefrac must be between 0 and 1"))

#! format: off
packmol_input_keywords = Dict{String,Function}(
    "filetype"                 => (T, val) -> (:filetype, _parse_value(String, "filetype", val)),
    "output"                   => (T, val) -> (:output_file, _parse_value(String, "output", val)),
    "tolerance"                => (T, val) -> (:tolerance, _parse_value(T, "tolerance", val)),
    "tolerance_precision"      => (T, val) -> (:tolerance_precision, _parse_value(T, "tolerance_precision", val)),
    "constraint_precision"     => (T, val) -> (:constraint_precision, _parse_value(T, "constraint_precision", val)),
    "maxit"                    => (T, val) -> (:maxit, _parse_value(Int, "max_iter", val)),
    "add_amber_ter"            => (T, val) -> (:add_amber_ter, _parse_value(Bool, "add_amber_ter", val)),
    "amber_ter_preserve"       => (T, val) -> (:amber_ter_preserve, true),
    "add_box_sides"            => (T, val) -> (:add_box_sides, true),
    "connect"                  => (T, val) -> (:connect, _parse_options(String, "connect", val, ("yes" => true, "no" => false))),
    "randominitialpoint"       => (T, val) -> (:randominitialpoint, true),
    "seed"                     => (T, val) -> (:seed, _parse_value(Int, "seed", val)),
    "avoid_overlap"            => (T, val) -> (:avoid_overlap, _parse_options(String, "avoid_overlap", val, ("yes" => true, "no" => false))),
    "writeout"                 => (T, val) -> (:writeout, _parse_value(Int, "writeout", val)),
    "writebad"                 => (T, val) -> (:writebad, true),
    "optimization_print_level" => (T, val) -> (:optimization_print_level, _parse_value(Int, "optim_print_level", val)),
    "chkgrad"                  => (T, val) -> (:check_gradient, true),
    "use_short_tol"            => (T, val) -> (:use_short_tol, true),
    "short_tol_dist"           => (T, val) -> (:short_tol_dist, _parse_value(T, "short_tol_dist", val)),
    "short_tol_scale"          => (T, val) -> (:short_tol_scale, _parse_value(T, "short_tol_scale", val)),
    "short_radius"             => (T, val) -> (:short_radiues, _parse_value(T, "short_radius", val)),
    "short_radius_scale"       => (T, val) -> (:short_radius_scale, _parse_value(T, "short_radius_scale", val)),
    "movebadrandom"            => (T, val) -> (:movebadrandom, true),
    "movefrac"                 => (T, val) -> (:movefrac, _parse_value(T, "movefrac", val; _val_check=_check_movefrac)),
)
#! format: on

packmol_legacy_keywords = Dict{String,String}(
    "fscale" => "fscale legacy keyword was ignored.",
    "fbins" => "fbins legacy keyword was ignored.",
    "iprint1" => "iprint1 legacy keyword was ignored, instead use: optim_print_level",
    "iprint2" => "iprint1 legacy keyword was ignored, instead use: optim_print_level",
    "iprint3" => "iprint1 legacy keyword was ignored, instead use: optim_print_level",
    "precision" => "precision legacy keyword was ignored, instead use: tolerance_precision and/or constraint_precision",
)

#=

    read_packmol_input

Reads the packmol input file to create a `PackmolSystem` object.

=#
function read_packmol_input(input_file::String; D::Int=3, T::DataType=Float64)
    input_data = Dict{Symbol,Any}(
        :input_file => input_file,
        :structure_types => StructureType{D,T}[],
        :atoms => AtomData{D,T}[]
        :molecule_positions => MoleculePosition{D,T}[]
        :nmols => 0,
    )
    structure_section = false
    open(input_file) do io
        for line in eachline(io)
            line = strip(line)
            (startswith(line, "#") || isempty(line)) && continue
            keyword, values... = split(line)
            if haskey(packmol_input_keywords, keyword)
                field_name, field_value = packmol_input_keywords[keyword](T, first(values))
                input_data[field_name] = field_value
                continue
            elseif haskey(packmol_legacy_keywords, keyword)
                @warn begin
                    packmol_legacy_keywords[keyword]
                end _line = line _file = input_file
                continue
            end
            if keyword == "structure"
                structure_section = true
                continue
            end
            if keyword == "end" && values[1] == "structure"
                if structure_section
                    structure_section = false
                end 
                continue
            end
            if structure_section
                continue
            end
            throw(ArgumentError("Keyword $keyword not recognized"))
        end
    end
    if structure_section
        throw(ArgumentError("Structure block not closed"))
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
    mol_index = 0
    atom_index = 0
    atoms = input_data[:atoms]
    molecule_positions = input_data[:molecule_positions]
    for (itype, structure_type) in enumerate(input_data[:structure_types])
        for _ in 1:structure_type.number_of_molecules
            mol_index += 1 
            push!(molecule_positions, 
                MoleculePosition(zeros(SVector{D,T}), zeros(SVectors{D,T}))
            )
            for iatom in 1:structure_type.number_of_atoms
                atom_index += 1
                push!(atoms, 
                    AtomData(
                        mol_index, 
                        itype,
                        structure_type.radii[iatom],
                        structure_type.constraints[iatom]
                    )
                )
            end
        end
    end
    input_data[:nmols] = mol_index

    return PackmolSystem{D,T}(; input_data...)
end

@testitem "_parse_value" begin
    using Packmol: _parse_value
    @test _parse_value(Int, "max_iter", "100") == 100
    @test _parse_value(Float64, "tolerance", "2.0") == 2.0
    @test _parse_value(Float32, "tolerance", "2.0") == 2.0f0
    @test_throws ArgumentError _parse_value(Int, "max_iter", "100.0")
end

@testitem "_parse_options" begin
    using Packmol: _parse_options
    @test _parse_options(String, "connect", "yes", ("yes" => true, "no" => false)) == true
    @test _parse_options(String, "connect", "no", ("yes" => true, "no" => false)) == false
    @test_throws ArgumentError _parse_options(String, "connect", "nop", ("yes" => true, "no" => false))
end

@testitem "read_packmol_input" begin
    using Packmol: read_packmol_input
    file = Packmol.src_dir * "/../test/run_packmol/water_box.inp"
    sys = read_packmol_input(file)


    sys = read_packmol_input(file; T=Float32)

    sys = read_packmol_input(file; D=2)

end