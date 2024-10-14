struct AtomData{T}
    molecule_index::Int
    structure_type_index::Int
    radius::T
    constraints::Vector{Int}
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
struct MoleculePosition{D,T}
    cm::SVector{D,T}
    angles::SVector{D,T}
end
Base.copy(x::MoleculePosition) = MoleculePosition(x.cm, x.angles)
import Base: + 
+(x::MoleculePosition, y::MoleculePosition) = MoleculePosition(x.cm + y.cm, x.angles + y.angles)
Base.zero(::Type{MoleculePosition{D,T}}) where {D,T} = MoleculePosition(zero(SVector{D,T}), zero(SVector{D,T}))
MoleculePosition(x,y,z,β,γ,θ) = MoleculePosition(SVector(x,y,z), SVector(β,γ,θ))

struct FixedMoleculeData{D,T}
    fixed::Bool
    position::MoleculePosition{D,T}
end
Base.zero(::Type{FixedMoleculeData{D,T}}) where {D,T} = FixedMoleculeData(false, zero(MoleculePosition{D,T}))
