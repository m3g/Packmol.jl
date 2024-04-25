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
# Dictionary containing the functions required for parsing each type
# of constraint. The functions for parsing each constraint are defined
# in the corresponding files, and appended to this dictionary.
#
const parse_constraint = Dict{String,Function}()



