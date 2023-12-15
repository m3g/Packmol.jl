module Packmol

using TestItems
using Parameters
using StaticArrays
import CellListMap
using CellListMap.PeriodicSystems
using SPGBox
import LinearAlgebra: norm
import Statistics: mean

const src_dir = @__DIR__

# API: exported functions
export pack_monoatomic!

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

# Monatomic packing
include("./mono_atomic.jl")

# Runner for the legacy packmol 
include("./packmol_runner.jl")

end
