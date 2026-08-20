#=

This functions apply the chain rule to compute the gradient of the minimum-distance
function relative to the rotations and translation of the rigid-body molecules, from 
the gradient computed from the displacement of the atoms in cartesian coordinates. 

This was implemented originally in the Fortran code of Packmol, by J. M. Martínez.

>> The function modifies the fg.g vector

=#
function chain_rule!(fg, packmol_system::PackmolSystem{D,T}) where {D,T}
    imol_offset = 0
    iat_offset = 0
    for structure_type in packmol_system.structure_types
        natoms_st = structure_type.natoms
        nmols_st = structure_type.number_of_molecules
        st_imol_offset = imol_offset
        st_iat_offset = iat_offset
        @sync for irange in chunks(1:nmols_st; n=Threads.nthreads())
            @spawn for i in irange
                imol = st_imol_offset + i
                ifmol = st_iat_offset + (i - 1) * natoms_st + 1
                ilmol = ifmol + natoms_st - 1
                gmol = partial_derivatives(
                    fg.g[imol],
                    packmol_system.molecule_positions[imol].angles,
                    structure_type.reference_coordinates,
                    @view(fg.gxcar[ifmol:ilmol]),
                )
                fg.g[imol] = gmol
            end
        end
        imol_offset += nmols_st
        iat_offset += nmols_st * natoms_st
    end
    return packmol_system
end

# Same as chain_rule! above, restricted to an explicit, possibly
# non-contiguous list of molecules. Used by adjust_constraints!'s active-set
# loop, alongside compute_positions_and_constraints_for_mols! — only
# molecules with a fresh fg.gxcar this pass need their gradient chain-ruled.
#
# partial_derivatives accumulates onto its `gmol` argument (see gcm/grot
# below), so chain_rule! normally relies on the caller having zeroed the
# *entire* fg.g array first (CellListMap.reset_output! in the full-system
# path). Here that would reintroduce an O(nmols) cost on every call, which is
# exactly what the active-set restriction is meant to avoid — so this passes
# a fresh zero seed per molecule directly instead of reading (possibly stale)
# fg.g[imol].
function chain_rule_for_mols!(
    fg, packmol_system::PackmolSystem{D,T},
    mol_list::Vector{Int}, mol_structure_type::Vector{Int}, mol_iat_first::Vector{Int},
) where {D,T}
    @sync for (_, mrange) in enumerate(chunks(1:length(mol_list); n=Threads.nthreads()))
        @spawn for k in mrange
            imol = mol_list[k]
            structure_type = packmol_system.structure_types[mol_structure_type[imol]]
            ifmol = mol_iat_first[imol]
            ilmol = ifmol + structure_type.natoms - 1
            gmol = partial_derivatives(
                zero(MoleculePosition{D,T}),
                packmol_system.molecule_positions[imol].angles,
                structure_type.reference_coordinates,
                @view(fg.gxcar[ifmol:ilmol]),
            )
            fg.g[imol] = gmol
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
        # dx_atom/d_angle_d = (∂R/∂angle_d) * x_ref
        # df/d_angle_d = sum_i dot(gx_car_i, (∂R/∂angle_d) * x_ref_i)
        grot += SVector(ntuple(d -> dot(gx, ∂v_matrices[d] * x), D))
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
    # Must match eulermat(β, γ, θ) from rigid_body.jl exactly.
    # eulermat uses: c1=cos(β), s1=sin(β), c2=cos(γ), s2=sin(γ), c3=cos(θ), s3=sin(θ)
    # Here: sb=s1, cb=c1, sg=s2, cg=c2, st=s3, ct=c3
    (sb, cb), (sg, cg), (st, ct) = sincos.(angles)
    #!format: off
    # ∂R/∂β: only rows 2,3 are nonzero (row 1 has no β dependence)
    ∂v∂β = @SMatrix[
             zero(T)            zero(T)   zero(T)
        -sb*st+cb*ct*sg   -sb*ct-cb*sg*st   -cb*cg
         cb*st+sb*ct*sg   -sb*sg*st+cb*ct   -sb*cg
    ]
    # ∂R/∂γ: all rows have γ dependence
    ∂v∂γ = @SMatrix[
        -sg*ct        sg*st      cg
         sb*cg*ct    -sb*cg*st   sb*sg
        -cb*cg*ct     cb*cg*st  -cb*sg
    ]
    # ∂R/∂θ: column 3 is zero (no θ dependence in column 3 of eulermat)
    ∂v∂θ = @SMatrix[
              -cg*st            -cg*ct    zero(T)
         cb*ct-sb*sg*st   -cb*st-sb*sg*ct zero(T)
         sb*ct+cb*sg*st    cb*sg*ct-sb*st  zero(T)
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
