#=

This functions apply the chain rule to compute the gradient of the minimum-distance
function relative to the rotations and translation of the rigid-body molecules, from 
the gradient computed from the displacement of the atoms in cartesian coordinates. 

This was implemented originally in the Fortran code of Packmol, by J. M. Martínez.

>> The function modifies the fg.g vector

=#
function chain_rule!(fg, packmol_system::PackmolSystem{D,T}) where {D,T}
    imol = 0
    iat = 0
    for structure_type in packmol_system.structure_types
        for _ in 1:structure_type.number_of_molecules
            imol += 1
            ifmol = iat + 1
            ilmol = iat + structure_type.natoms
            gmol = partial_derivatives(
                fg.g[imol],
                packmol_system.molecule_positon[imol].angles,
                structure_type.reference_coordinates, 
                @view(fg.gxcar[ifmol:ilmol]),
            )
            fg.g[imol] = gmol
            iat += structure_type.natoms
        end
    end
    return system
end

#=

For a single molecule, computes the partial derivatives of the gradient of the
minimum distance between atoms relative to the rotations and translations of the
molecule, given the gradient of the minimum distance between atoms relative to the
displacement of the atoms in cartesian coordinates.

This function does not mutate any argument. It returns the gradient for the 
molecule as a MoleculePosition data structure, to be used in the chain_rule! function.

=#
function partial_derivatives(
    gmol::MoleculePosition{D,T}, 
    angles::SVector{D,T},
    reference_coordinates::Vector{SVector{D,T}},
    gxcar::SVector{D,T},
) where {D,T}
    gcm = gmol.cm
    grot = gmol.angles
    (sb, cb), (sg, cg), (st, ct) = sincos.(angles)
    #!format: off
    ∂v∂β = sum(SMatrix[
        -cb*sg*ct-sb*cg  -cb*cg*ct+sb*sg     cb*st
        -sb*sg*ct+cb*cg  -sb*cg*ct-cb*sg     sb*st
                zero(T)          zero(T)   zero(T)
    ]; dims=1)
    ∂v∂γ = sum(@SMatrix[
        -sb*cg*ct-cb*sg   sb*sg*ct-cb*cg  zero(T)
         cb*cg*ct-sb*sg  -sg*cb*ct-cg*sb  zero(T)
                  cg*st           -sg*st  zero(T)
    ]; dims=1)
    ∂v∂θ = sum(@SMatrix[
       -sb*sg*st  -sb*cg*st   -sb*ct
        cb*sg*st   cb*cg*st    cb*ct
           sg*ct      cg*ct      -st
    ]; dims=1)
    ∂rot = vcat(∂v∂β, ∂v∂γ, ∂v∂θ)
    #!format: on
    for i in eachindex(reference_coordinates, gxcar)
        x = reference_coordinates[i]
        gx = gxcar[i]
        gcm += gx
        grot += x' * ∂rot * gx
    end
    return MoleculePosition{D,T}(gcm, grot)
end

#=

2D case

=#
function partial_derivatives!(g, x::SVector{2}, system::PackmolSystem{2})
    error("not implemented.")
end

@testitem "gradient chain rule 3D" begin

end

@testitem "gradient chain rule 2D" begin

end
