#=

This functions apply the chain rule to compute the gradient of the minimum-distance
function relative to the rotations and translation of the rigid-body molecules, from 
the gradient computed from the displacement of the atoms in cartesian coordinates. 

This was implemented originally in the Fortran code of Packmol, by J. M. Martínez.

Extended here for the 2D case.

=#
function chain_rule!(packmol_system::PackmolSystem{D,T}, fg) where {D,T}
    fill!(zero(MoleculePositions{D,T}), packmol_system.gradient)
    imol = 0
    iat = 0
    for structure_type in packmol_system.structure_types
        for _ in 1:structure_type.number_of_molecules
            imol += 1
            ifmol = iat + 1
            ilmol = iat + structure_type.natoms
            partial_derivatives!(imol, structure_type, packmol_system, @view(fg.gxcar[ifmol:ilmol]))
            iat += structure_type.natoms
        end
    end
    return system
end

#=

3D case

=#
function partial_derivatives!(packmol_system, imol::Int, structure_type, gxcar)
    g_cm = packmol_system.gradient[imol].cm # vector of gradient relative to the center of mass
    g_rot = packmol_system.gradient[imol].angles # vector of gradient relative to the rotations
    (sb, cb), (sg, cg), (st, ct) = sincos.(packmol_system.molecule_position[imol].angles)
    #!format off
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
    #!format on
    for i in eachindex(structure_type.reference_coordinates, gxcar)
        x = structure_type.reference_coordinates[i]
        gx = gxcar[i]
        g_cm += gx
        g_rot += x' * ∂rot * gx
        packmol_system.gradient = MoleculePositions{D,T}(g_cm, g_rot)
    end
    return packmol_system.gradient
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
