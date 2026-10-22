#
# Types of constraints 
#
abstract type Constraint{Placement,N,T} end

struct Inside end
struct Outside end
struct Over end
struct Below end
export Inside, Outside, Over, Below

#
# Generic show function for constraints
#
function Base.show(io::IO, ::MIME"text/plain", c::Constraint)
    println(io, typeof(c))
    fnames = fieldnames(typeof(c))
    for i in 1:length(fnames)-1
        println(io, "    $(fnames[i]) = $(getfield(c,fnames[i]))")
    end
    print(io, "    $(fnames[end]) = $(getfield(c,fnames[end]))")
end

#
# Default weights
const weight_default = Dict{Symbol,Float64}()

#
# Parses the data of a constraint line of the input file, returning the
# constraint. Each input file keyword (e.g. "inside box") has its own method,
# `parse_constraint(::Val{Symbol(keyword)}, data, ::Type{T}, filename)`, defined
# in the corresponding file (see also `constraint_keywords` in constraint_types.jl).
#
function parse_constraint end



