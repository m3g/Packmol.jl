#
# The abstract type of the constraints, `Constraint{Placement,N,T}`, and the markers of
# their placement (`Inside`, `Outside`, `Over`, `Below`); the generic `show`, the default
# weights of the penalties (set by each constraint file), the generic `parse_constraint`,
# and the fallback of `_inside_box`.
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
# Default weights of the penalties, by constraint (set in the file of each constraint)
const weight_default = Dict{Symbol,Float64}()

#
# Parses the data of a constraint line of the input file, returning the
# constraint. Each input file keyword (e.g. "inside box") has its own method,
# `parse_constraint(::Val{Symbol(keyword)}, data, ::Type{T}, filename)`, defined
# in the corresponding file (see also `constraint_keywords` in constraint_types.jl).
#
function parse_constraint end

#
# The bounding box `(lo, hi)` of the region of an "inside" constraint, defined for the
# constraints with a bounded region in their files; `nothing` for the others. It bounds
# the sampling of the initial placement, the Monte Carlo points of `free_volume`, and the
# grid of `voxel_place`.
#
_inside_box(::Constraint) = nothing
