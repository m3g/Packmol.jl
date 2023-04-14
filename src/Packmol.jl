module Packmol

using TestItems
using Parameters
using StaticArrays
using CellListMap
using SPGBox
using LinearAlgebra: norm

# Flag for internal function doc entries
const INTERNAL = "Internal function or structure - interface may change."

include("./rigid_body.jl")
include("./constraints.jl")

@with_kw struct Structure{N,T}
    number::Int
    coor::Vector{SVector{N,T}}
    fixed::Bool = false
    radii::Vector{T} = [ T(2) for _ in 1:length(coor) ]
    residue_numbering::Int = 1
    connect::Vector{Vector{Int}} = Vector{Int}[]
    constraints::Vector{<:Constraint}
end

end
