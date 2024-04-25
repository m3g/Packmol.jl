module Packmol

using TestItems: @testitem
using StaticArrays: SVector, SMatrix, @SMatrix
using LinearAlgebra: norm
using Statistics: mean
using Base: @kwdef
using Base.Threads: @spawn

const src_dir = @__DIR__

# API: exported functions
export pack_monoatomic!

# Rigid body transformations 
include("./rigid_body/rigid_body.jl")
include("./rigid_body/chain_rule.jl")

# Constraints
include("./constraints/constraints_base.jl")
include("./constraints/boxes.jl")
include("./constraints/spheres.jl")

@kwdef struct Structure{N,T}
    number::Int
    coor::Vector{SVector{N,T}}
    fixed::Bool = false
    radii::Vector{T} = [ T(2) for _ in 1:length(coor) ]
    residue_numbering::Int = 1
    connect::Vector{Vector{Int}} = Vector{Int}[]
    constraints::Vector{<:Constraint}
end

@kwdef PackmolSystem{N,T}

end

# Monatomic packing
include("./mono_atomic.jl")

# Runner for the legacy packmol 
include("./packmol_runner.jl")

end
