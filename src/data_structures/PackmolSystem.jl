@kwdef mutable struct PackmolSystem{D,T}
    input_file::String
    output_file::String = "packmol_output.pdb"
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
end

function Base.show(io::IO, sys::PackmolSystem{D,T}) where {D,T}
    println(io, "PackmolSystem{D,T} with $(length(sys.structure_types)) structure types")
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
    value = tryparse(T, input_value) 
    if isnothing(value)
        throw(ArgumentError("""\n
            Could not parse value $input_value for keyword $keyword. Expected type: $T.

        """))
    end
    _val_check(value)
    return value
end

function _parse_options(T::DataType, keyword::String, input_value, options::Tuple) 
    input_value = _parse_value(T, keyword, input_value)
    for (option, value) in options
        input_value == option && return value
    end
    throw(ArgumentError("""\n 
        Could not parse value "$input_value" for keyword: $keyword. 
        Expected one of: $(join(options, ','))

    """))
end

_check_movefrac(x) = 0.0 <= x <= 1.0 ? x : throw(ArgumentError("movefrac must be between 0 and 1"))

##! format: off
packmol_input_keywords = Dict{String, Function}(
    "output"                   => (T, val) -> _parse_value(String, "output", val),
    "tolerance"                => (T, val) -> _parse_value(T, "tolerance", val),
    "tolerance_precision"      => (T, val) -> _parse_value(T, "tolerance_precision", val),
    "constraint_precision"     => (T, val) -> _parse_value(T, "constraint_precision", val),
#    "precision"                => (Real,   [(:tolerance_precision => :first), (:constraint_precision => :first)]),
    "maxit"                 => (T, val) -> _parse_value(Int, "max_iter", val),
    "add_amber_ter"            => (T, val) -> _parse_value(Bool, "add_amber_ter", val),
    "amber_ter_preserve"       => (T, val) -> true,
    "add_box_sides"            => (T, val) -> true,
    "connect"                  => (T, val) -> _parse_options(String, "connect", val, ("yes" => true, "no" => false)),
    "randominitialpoint"       => (T, val) -> true,
    "seed"                     => (T, val) -> _parse_value(Int, "seed", val),
    "avoid_overlap"            => (T, val) -> _parse_options(String, "avoid_overlap", val, ("yes" => true, "no" => false)),
    "writeout"                 => (T, val) -> _parse_value(Int, "writeout", val),
    "writebad"                 => (T, val) -> true,
    "optimization_print_level" => (T, val) -> _parse_value(Int, "optim_print_level", val),
    "chkgrad"                  => (T, val) -> true,
    "use_short_tol"            => (T, val) -> true,
    "short_tol_dist"           => (T, val) -> _parse_value(T, "short_tol_dist", val),
    "short_tol_scale"          => (T, val) -> _parse_value(T, "short_tol_scale", val),
    "short_radius"             => (T, val) -> _parse_value(T, "short_radius", val),
    "short_radius_scale"       => (T, val) -> _parse_value(T, "short_radius_scale", val),
    "movebadrandom"            => (T, val) -> true,
    "movefrac"                 => (T, val) -> _parse_value(T, "movefrac", val; _val_check=_check_movefrac),
)
#! format: on

packmol_legacy_keywords = Dict{String, String}(
    "fscale" => "fscale legacy keyword was ignored.",
    "fbins" => "fbins legacy keyword was ignored.",
    "iprint1" => "iprint1 legacy keyword was ignored, instead use: optim_print_level",
    "iprint2" => "iprint1 legacy keyword was ignored, instead use: optim_print_level",
    "iprint3" => "iprint1 legacy keyword was ignored, instead use: optim_print_level",
)

function set_packmol_input_data!(input_data, line; keywords = packmol_input_keywords)
    keyword, values = split(line)
    if haskey(keyword, keywords)
        for key in keywords[keyword]
            input_data[key] = values[1]
        end
    else
        throw(ArgumentError("Keyword $keyword not recognized"))
    end
    return input_data
end



#=

    read_packmol_input

Reads the packmol input file to create a `PackmolSystem` object.

=#
function read_packmol_input(input_file::String; D=3, T=Float64)
    packmol_input_data = Dict{Symbol,Any}(
        :input_file => input_file,
    )
    open(input_file) do io
        structure_field = 0
        for line in eachline(io)
            line = strip(line)
            (startswith(line, "#") || isempty(line)) && continue
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
    end
    return PackmolSystem{D,T}(;packmol_input_data...)
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
    file = Packmol.src_dir*"/../test/run_packmol/water_box.inp"
    sys = read_packmol_input(file)


    sys = read_packmol_input(file; T = Float32)

    sys = read_packmol_input(file; D=2)

end
