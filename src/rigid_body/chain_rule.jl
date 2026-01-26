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
                packmol_system.molecule_positions[imol].angles,
                structure_type.reference_coordinates, 
                @view(fg.gxcar[ifmol:ilmol]),
            )
            fg.g[imol] = gmol
            iat += structure_type.natoms
        end
    end
    return packmol_system
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
    gxcar::AbstractVector{SVector{D,T}},
) where {D,T}
    # Compute rotation derivative matrices
    ∂v_matrices = rotation_derivative_matrices(angles)
    
    # Initialize gradients with current values
    gcm = gmol.cm
    grot = gmol.angles
    
    # Loop over all atoms in the molecule
    for i in eachindex(reference_coordinates, gxcar)
        x = reference_coordinates[i]
        gx = gxcar[i]
        # Gradient contribution from center of mass translation
        gcm += gx
        # Gradient contribution from rotations (chain rule)
        # For each angle, compute: sum_k (dx/dangle)_k * gx_k
        # where (dx/dangle)_k = sum_j x_j * (dv_j/dangle)_k
        grot += SVector(ntuple(d -> dot(x, ∂v_matrices[d] * gx), D))
    end
    return MoleculePosition{D,T}(gcm, grot)
end

#=

Compute the derivative matrices of the rotation matrix with respect to each angle.
Returns a tuple of D matrices, where D is the dimension (2 or 3).

For 2D: Returns (∂R/∂θ,) where θ is the single rotation angle
For 3D: Returns (∂R/∂β, ∂R/∂γ, ∂R/∂θ) for the three Euler angles

=#
function rotation_derivative_matrices(angles::SVector{3,T}) where {T}
    # 3D case: Euler angles (β, γ, θ)
    (sb, cb), (sg, cg), (st, ct) = sincos.(angles)
    #!format: off
    # Partial derivatives of the rotation matrix columns with respect to each Euler angle
    # These are 3x3 matrices where each column is dv_i/dangle
    ∂v∂β = @SMatrix[
        -cb*sg*ct-sb*cg  -sb*sg*ct+cb*cg  cb*st
        -cb*cg*ct+sb*sg  -sb*cg*ct-cb*sg  sb*st
                zero(T)          zero(T)  zero(T)
    ]
    ∂v∂γ = @SMatrix[
        -sb*cg*ct-cb*sg   cb*cg*ct-sb*sg  zero(T)
         sb*sg*ct-cb*cg  -sg*cb*ct-cg*sb  zero(T)
                  cg*st           -sg*st  zero(T)
    ]
    ∂v∂θ = @SMatrix[
        sb*sg*st   cb*sg*st   sg*ct
        sb*cg*st   cb*cg*st   cg*ct
       -sb*ct     -cb*ct     -st
    ]
    #!format: on
    return (∂v∂β, ∂v∂γ, ∂v∂θ)
end

function rotation_derivative_matrices(angles::SVector{2,T}) where {T}
    # 2D case: Single rotation angle θ
    θ = angles[1]
    sθ, cθ = sincos(θ)
    #!format: off
    # Derivative of 2D rotation matrix with respect to θ
    # R(θ) = [cos(θ) -sin(θ)]    =>    dR/dθ = [-sin(θ) -cos(θ)]
    #        [sin(θ)  cos(θ)]                   [ cos(θ) -sin(θ)]
    ∂v∂θ = @SMatrix[
        -sθ  cθ
        -cθ -sθ
    ]
    #!format: on
    return (∂v∂θ,)
end

@testitem "gradient chain rule 3D" begin

end

@testitem "gradient chain rule 2D" begin

end
